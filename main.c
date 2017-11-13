#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <termios.h>
#include <stdlib.h>
#include <pthread.h>

#include "FileStruct.h"
#include "Queue.h"
#include "NewKalmanFilter.h"
#include "MahonyAHRS.h"
#include "parameter.h"
#include "Turnning.h"

#define accCalNum 500
#define windowWidth 12
#define rad2deg 57.29578

struct Quaternion euarToQuaternion(double yaw, double pitch, double roll);

void quaternionToMahonyR(struct Quaternion Q);

void accToFn(double x, double y, double z);

void *start_deal(void *que);

void write_raw_data(struct Data d, FILE *fp);

void write_slidefit_data(struct slideRes r, FILE *fp);

void write_raw_gps(struct Data d, FILE *fp);

void setMahonyR();


struct STime stcTime;
struct SAcc stcAcc;
struct SGyro stcGyro;
struct SAngle stcAngle;
struct SMag stcMag;
struct SLonLat stcLonLat;
struct SGPSV stcGPSV;
struct SQUAT stcQUAT;
struct SDOP stcDOP;

//Re长半轴 r短半轴 f椭-球扁率 e椭球偏心率 wie地球自转角速率
const double earthRe = 6378137, earthr = 6356752.3142, earthf = 1 / 298.257, earthe = 0.0818, earthwie = 7.292e-5;
struct Data input_data;
struct CircleQueue queue;
double accCalStc_X = -0.09092, accCalStc_Y = 0.081208, accCalStc_Z = 0.015632;
double pitch0 = 0, roll0 = 0, yaw0 = 0;
double Pitch = 0, Roll = 0, Yaw = 0;
struct Quaternion Q0;
extern float q0, q1, q2, q3;

int count_itr = 0;
const int worksize = 1500;
//-------------------deal thread set-------------------------------------------
int datacnt = 0;
double samplePeriod = 0.1;
const double INS_drift_weight = 0;
int doTurning = 1;//是否用turning算法航位推算覆盖Vccq,0为否
double stepP[3] = {0,0,0};
double gyoP[3] = {0,0,0};
double magP[3] = {0,0,0};


int firstGPSOff = 0;
int GPSOff = 0;
//定义全局的旋转矩阵
double *mahonyR;
double m_temp[3];
const int width = 12;//平滑窗口大小
int firstSmooth = 0;

double Px = 0, Py = 0;

double *Fn;
double Vccq[3];
double Vcc[3];
double lastVx = 0.0, lastVy = 0.0, lastVz = 0.0;

double Rm = 0, Rn = 0, R0 = 0;
double Last_L = 0.0, Last_E = 0.0, Last_h = 0.0;
double pre_L = 0.0, pre_E = 0.0, pre_h = 0.0;
double INS_L = 0.0, INS_E = 0.0, INS_h = 0.0;
double speedE = 0, speedN = 0, speedH = 0;
extern double G0;
double tao = 0;
struct Position AddS = {0, 0, 0, 0};
double*** smoothRes;
struct slideRes sum_buf;

void slideFilter(struct Data d) {
    if (firstSmooth == 0) {
        for (int i = 0; i < width; i++) {
            smoothRes[0][0][i] = d.acc_x;
            smoothRes[0][1][i] = d.acc_y;
            smoothRes[0][2][i] = d.acc_z;
            smoothRes[1][0][i] = d.gyo_x;
            smoothRes[1][1][i] = d.gyo_y;
            smoothRes[1][2][i] = d.gyo_z;
            smoothRes[2][0][i] = d.m_x;
            smoothRes[2][1][i] = d.m_y;
            smoothRes[2][2][i] = d.m_z;
            smoothRes[3][0][i] = d.Yaw;
            smoothRes[3][1][i] = d.GPSV;
        }
        firstSmooth = 1;
    } else {
        for (int i = 1; i < width; i++) {
            smoothRes[0][0][i - 1] = smoothRes[0][0][i];
            smoothRes[0][1][i - 1] = smoothRes[0][1][i];
            smoothRes[0][2][i - 1] = smoothRes[0][2][i];
            smoothRes[1][0][i - 1] = smoothRes[1][0][i];
            smoothRes[1][1][i - 1] = smoothRes[1][1][i];
            smoothRes[1][2][i - 1] = smoothRes[1][2][i];
            smoothRes[2][0][i - 1] = smoothRes[2][0][i];
            smoothRes[2][1][i - 1] = smoothRes[2][1][i];
            smoothRes[2][2][i - 1] = smoothRes[2][2][i];
            smoothRes[3][0][i - 1] = smoothRes[3][0][i];
            smoothRes[3][1][i - 1] = smoothRes[3][1][i];
        }
        smoothRes[0][0][width - 1] = d.acc_x;
        smoothRes[0][1][width - 1] = d.acc_y;
        smoothRes[0][2][width - 1] = d.acc_z;
        smoothRes[1][0][width - 1] = d.gyo_x;
        smoothRes[1][1][width - 1] = d.gyo_y;
        smoothRes[1][2][width - 1] = d.gyo_z;
        smoothRes[2][0][width - 1] = d.m_x;
        smoothRes[2][1][width - 1] = d.m_y;
        smoothRes[2][2][width - 1] = d.m_z;
        smoothRes[3][0][width - 1] = d.Yaw;;
        smoothRes[3][1][width - 1] = d.GPSV;
    }
    for(int j = 0;j < width;j++){
        sum_buf.slideAcc_x += smoothRes[0][0][j];
        sum_buf.slideAcc_y += smoothRes[0][1][j];
        sum_buf.slideAcc_z += smoothRes[0][2][j];
        sum_buf.slideGyo_x += smoothRes[1][0][j];
        sum_buf.slideGyo_y += smoothRes[1][1][j];
        sum_buf.slideGyo_z += smoothRes[1][2][j];
        sum_buf.slideMag_x += smoothRes[2][0][j];
        sum_buf.slideMag_y += smoothRes[2][1][j];
        sum_buf.slideMag_z += smoothRes[2][2][j];
        sum_buf.slideGps_yaw += smoothRes[3][0][j];
        sum_buf.slideGps_v += smoothRes[3][1][j];
    }
    sum_buf.slideAcc_x /= width;
    sum_buf.slideAcc_y /= width;
    sum_buf.slideAcc_z /= width;
    sum_buf.slideGyo_x /= width;
    sum_buf.slideGyo_y /= width;
    sum_buf.slideGyo_z /= width;
    sum_buf.slideMag_x /= width;
    sum_buf.slideMag_y /= width;
    sum_buf.slideMag_z /= width;
    sum_buf.slideGps_yaw /= width;
    sum_buf.slideGps_v /= width;
}

void *start_deal(void *que) {
    int firstQuatCalStatus = 1;
    double accD0_x = 0, accD0_y = 0, accD0_z = 0;
    double magD0_x = 0, magD0_y = 0, magD0_z = 0;
    struct CircleQueue *queue = (struct CircleQueue *) que;
    FILE *f_raw_data, *f_slide, *f_gps, *f_imu_a, *f_imu_v, *f_imu_p,*f_tokf,*f_kfresult,*fresult;

    f_raw_data = fopen("raw_data.csv", "w");
    f_slide = fopen("smoothdata.csv", "w");
    f_gps = fopen("gps.csv", "w");
    f_imu_a = fopen("imu_fn.csv", "w");
    f_imu_v = fopen("imu_velocity_ENU.csv", "w");
    f_imu_p = fopen("imu_position_ENU.csv", "w");
    fresult = fopen("result.csv", "w");
    // f_tokf = fopen("file_tokf.csv","w");
    // f_kfresult = fopen("file_kfresult.csv","w");

    while (1) {
        printf("in the deal thread!\n");
        while (isEmpty(queue)) {}
        struct Data data = DeQueue(queue);

        double ax = -(data.acc_x - accCalStc_X);
        double ay = -(data.acc_y - accCalStc_Y);
        double az = -(data.acc_z - accCalStc_Z);
        double gx = data.gyo_x;
        double gy = data.gyo_y;
        double gz = data.gyo_z;
        double mx = data.m_x;
        double my = data.m_y;
        double mz = -data.m_z;


        if (firstGPSOff == 0) {
            if (data.SN < 4) {
                printf("Please wait until GPS is working..\r\n%f\r\n\r\n", accCalStc_Z);
                continue;
            } else if (data.Lattitude != 0) {
                //第一次接收到GPS信号，确定起点经纬度，或者是GPS信号丢失后重新获取，然后重新设置起点
                write_raw_gps(data, f_gps);

                Last_E = data.Longitude / rad2deg;
                Last_L = data.Lattitude / rad2deg;
                Last_h = data.Height;
                INS_E = Last_E;
                INS_L = Last_L;
                INS_h = Last_h;

                pre_E = Last_E;
                pre_L = Last_L;
                pre_h = Last_h;
                printf("Now we set the GPS first location: E:%.5f L:%.5f...............session 0\n", Last_E * rad2deg,
                       Last_L * rad2deg);
                firstGPSOff = 1;
                datacnt++;
                
				accD0_x = ax;
				accD0_y = ay;
				accD0_z = az;
				magD0_x = mx;
				magD0_y = my;
				magD0_z = mz;

				pitch0 = atan2(-accD0_y, -accD0_z);
				roll0 = atan2(accD0_x, -accD0_z);
				yaw0 = atan2(-magD0_y * cos(roll0) + magD0_z * sin(roll0),
							 magD0_x * cos(pitch0) + magD0_y * sin(pitch0) * sin(roll0) -
							 magD0_z * sin(pitch0) * cos(roll0));

				pitch0 = pitch0 * rad2deg;
				roll0 = roll0 * rad2deg;
				yaw0 = -yaw0 * rad2deg;

				Q0 = euarToQuaternion(yaw0, pitch0, roll0);
				//将第一组四元数传入到mahony更新方法里并设置好初始的旋转矩阵
				quaternionToMahonyR(Q0);
                   
               
            }
        } else if (firstGPSOff == 1) {
            write_raw_data(data, f_raw_data);

            slideFilter(data);
            write_slidefit_data(sum_buf, f_slide);

            speedN = sum_buf.slideGps_v * cos(sum_buf.slideGps_yaw * 3.1415926 / 180) / 3.6;
            speedE = sum_buf.slideGps_v * sin(sum_buf.slideGps_yaw * 3.1415926 / 180) / 3.6;
            speedH = data.Height - Last_h;//GPS速度
            //更新四元数
            MahonyAHRSupdate((float) (gx / rad2deg), (float) (gy / rad2deg), (float) (gz / rad2deg),
                             (float)ax, (float)ay, (float)az,
                             (float)mx, (float)my, (float)mz);
            Yaw = -atan2(2 * q1 * q2 + 2 * q0 * q3, 2 * q0 * q0 + 2 * q0 * q2 - 1) * rad2deg;
            Pitch = asin(2 * q2 * q3 + 2 * q0 * q1) * rad2deg;
            Roll = -atan2(2 * q1 * q3 + 2 * q0 * q2, 2 * q0 * q0 + 2 * q3 * q3 - 1) * rad2deg;
            //更新旋转矩阵
            setMahonyR();

            double *resultOrientation = TurnningTest(gx, gy, gz, -mx, -my, mz);
            stepP[0] += sin(resultOrientation[3] * 3.1415926 / 180);
            stepP[1] += cos(resultOrientation[3] * 3.1415926 / 180);
            gyoP[0] += sin(resultOrientation[0] * 3.1415926 / 180);
            gyoP[1] += cos(resultOrientation[0] * 3.1415926 / 180);
            magP[0] += sin(resultOrientation[2] * 3.1415926 / 180);
            magP[1] += cos(resultOrientation[2] * 3.1415926 / 180);

            //测试单步长匀速路径
            Px += cos(Yaw * 3.1415926 / 180);
            Py += sin(Yaw * 3.1415926 / 180);

            //计算更新的子午曲率半径Rm和卯酉曲率半径Rn以及曲率平均半径R0
            Rm = earthRe * (1 - 2 * earthf + 3 * earthf * sin(Last_L) * sin(Last_L));
            Rn = earthRe * (1 + earthf * sin(Last_L) * sin(Last_L));
            R0 = sqrt(Rm * Rm + Rn * Rn);

            //旋转加速度到东北天导航坐标系进而计算绝对速度微分Vccq
            accToFn(data.acc_x - accCalStc_X, data.acc_y - accCalStc_Y, data.acc_z - accCalStc_Z);
            Fn = MatMulk(Fn, 3, 1, G0);
            Vccq[0] = -Fn[0];
            Vccq[1] = Fn[1];
            Vccq[2] = Fn[2] - G0;
            if(doTurning == 1){
                Vccq[0] = ay * G0 * sin(resultOrientation[3] * 3.1415926 / 180);
                Vccq[1] = ay * G0 * cos(resultOrientation[3] * 3.1415926 / 180);
            }

            Vcc[0] += Vccq[0] * samplePeriod;
            Vcc[1] += Vccq[1] * samplePeriod;
            Vcc[2] += Vccq[2] * samplePeriod;

            INS_L = INS_L + (lastVy / (Rm + Last_h)) * samplePeriod;
            INS_E = INS_E + (lastVx / (cos(Last_L) * (Rn + Last_h))) * samplePeriod;
            INS_h = INS_h - lastVz * samplePeriod;


            lastVx = Vcc[0];
            lastVy = Vcc[1];
            lastVz = Vcc[2];
            Last_L = INS_L;
            Last_h = INS_h;

            pre_E = data.Longitude;
            pre_L = data.Lattitude;
            pre_h = data.Height;

            if (data.SN < 4) {
                printf("lost GPS..\nNow We are in INS Mode...................................session 3\n");
                GPSOff = 1;

                INS_L = INS_L - (lastVy / (Rm + Last_h)) * samplePeriod * INS_drift_weight;
                INS_E = INS_E - (lastVy / (cos(Last_L) * (Rn + Last_h))) * samplePeriod * INS_drift_weight;
                INS_h = INS_h - lastVz * samplePeriod * INS_drift_weight;
                AddS.x = INS_E;
                AddS.y = INS_L;
                AddS.z = INS_h;
            } else if (data.SN >= 4) {
                if(count_itr++ == worksize){
                    firstGPSOff = 0;
                    tao = 0;
                    setNULL();
                    count_itr = 0;
                    continue;
                }
                //gps/ins mode
                write_raw_gps(data, f_gps);
                printf("Now we are in the GPS/INS mode...........................................session 2\n");


                //初速度设置

                if (GPSOff == 1) {
                    INS_E = data.Longitude / rad2deg;
                    INS_L = data.Lattitude / rad2deg;
                    INS_h = data.Height;
                    printf("Now we are setting the GPS location again after GPS lost...........................................session 2_1\n");
                    GPSOff = 0;
                }
//            initnum(data.Lattitude, data.Height, speedN, speedE);
                tao += samplePeriod;
                double Dpv[6] = {INS_L * rad2deg - pre_L, INS_E * rad2deg - pre_E, INS_h - pre_h, Vcc[0] - speedE,
                                 Vcc[1] - speedN, Vcc[2] - speedH};

//                char datas_tokf[10000];
//                sprintf(datas_tokf,"%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
//                        datacnt,Dpv[0],Dpv[1],Dpv[2]
//                        ,Dpv[3],Dpv[4],Dpv[5]
//                        ,Vcc[0],Vcc[1],Vcc[2]
//                        ,Last_L,Last_h,Fn[0]
//                        ,Fn[1],Fn[2]
//                        ,tao,Rm,Rn
//                        ,mahonyR[0],mahonyR[1],mahonyR[2]
//                        ,mahonyR[3],mahonyR[4],mahonyR[5]
//                        ,mahonyR[6],mahonyR[7],mahonyR[8]);
//                fprintf(f_tokf,datas_tokf,10000);
                double *XX = kalman_GPS_INS_pv(Dpv, Vcc[0], Vcc[1], Vcc[2], Last_L, Last_h, mahonyR, Fn, tao, Rm, Rn);

//                char datas_kfresult[10000];
//                sprintf(datas_kfresult,"%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
//                        datacnt,XX[0],XX[1],XX[2],XX[3],XX[4],XX[5],XX[6],XX[7],XX[8],XX[9]
//                        ,XX[10],XX[11],XX[12],XX[13],XX[15]);
//                fprintf(f_kfresult,datas_kfresult,10000);

                Vcc[0] -= XX[3];
                Vcc[1] -= XX[4];
                Vcc[2] -= XX[5];
                lastVx = Vcc[0];
                lastVy = Vcc[1];
                lastVz = Vcc[2];
                INS_L -= 0.29 * XX[6];
                INS_E -= 0.32 * XX[7];
                INS_h -= XX[8];

                AddS.x = INS_E;
                AddS.y = INS_L;
                AddS.z = INS_h;
            }
            fprintf(fresult, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
                    datacnt, Roll, Pitch, Yaw, Vccq[0], Vccq[1], Vccq[2], Vcc[0], Vcc[1], Vcc[2], INS_E * rad2deg,
                    INS_L * rad2deg, INS_h,resultOrientation[3]);
            fprintf(f_imu_a, "%d,%f,%f,%f\n", datacnt, Vccq[0], Vccq[1], Vccq[2]);
            fprintf(f_imu_v, "%d,%f,%f,%f\n", datacnt, Vcc[0], Vcc[1], Vcc[2]);
            fprintf(f_imu_p, "%d,%f,%f,%f\n", datacnt, INS_E * rad2deg, INS_L * rad2deg, INS_h);
            datacnt++;
			sum_buf.slideAcc_x = 0;
			sum_buf.slideAcc_y = 0;
			sum_buf.slideAcc_z = 0;
			sum_buf.slideGyo_x = 0;
			sum_buf.slideGyo_y = 0;
			sum_buf.slideGyo_z = 0;
			sum_buf.slideMag_x = 0;
			sum_buf.slideMag_y = 0;
			sum_buf.slideMag_z = 0;
			sum_buf.slideGps_yaw = 0;
			sum_buf.slideGps_v = 0;
        }
    }
    fclose(fresult);
}

void accToFn(double x, double y, double z) {
    m_temp[0] = x;
    m_temp[1] = y;
    m_temp[2] = z;

    Fn = MatMul(mahonyR, 3, 3, m_temp, 3, 1);
}


void setMahonyR() {
    extern float q0, q1, q2, q3;

    double q0q0 = q0 * q0;
    double q1q1 = q1 * q1;
    double q2q2 = q2 * q2;
    double q3q3 = q3 * q3;
    double q0q1 = q0 * q1;
    double q0q2 = q0 * q2;
    double q0q3 = q0 * q3;
    double q1q2 = q1 * q2;
    double q1q3 = q1 * q3;
    double q2q3 = q2 * q3;
    //定义Cnb
    mahonyR[0] = q0q0 + q1q1 - q2q2 - q3q3;
    mahonyR[1] = 2 * (q1q2 - q0q3);
    mahonyR[2] = 2 * (q1q3 + q0q2);
    mahonyR[3] = 2 * (q1q2 + q0q3);
    mahonyR[4] = q0q0 - q1q1 + q2q2 - q3q3;
    mahonyR[5] = 2 * (q2q3 - q0q1);
    mahonyR[6] = 2 * (q1q3 - q0q2);
    mahonyR[7] = 2 * (q2q3 + q0q1);
    mahonyR[8] = q0q0 - q1q1 - q2q2 + q3q3;
}

void Eular2R(struct Data d) {
    double x = d.angle_x;
    double y = d.angle_y;
    double z = d.angle_z;

    mahonyR[0] = cos(y) * cos(z);
    mahonyR[1] = cos(z) * sin(x) * sin(y) - cos(x) * sin(z);
    mahonyR[2] = sin(x) * sin(z) + cos(x) * cos(z) * sin(y);
    mahonyR[3] = cos(y) * sin(z);
    mahonyR[4] = cos(x) * cos(z) + sin(x) * sin(y) * sin(z);
    mahonyR[5] = cos(x) * sin(y) * sin(z) - cos(z) * sin(x);
    mahonyR[6] = -sin(y);
    mahonyR[7] = cos(y) * sin(x);
    mahonyR[8] = cos(x) * cos(y);
}

void write_raw_data(struct Data d, FILE *fp) {
    fprintf(fp, "%f,%f,%f,%f,%f,%f,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d\n"
            , d.acc_x, d.acc_y, d.acc_z
            , d.gyo_x, d.gyo_y, d.gyo_z
            , d.m_x, d.m_y, d.m_z
            , d.angle_x, d.angle_y, d.angle_z,
            d.Longitude, d.Lattitude, d.Height
            , d.Yaw, d.GPSV, d.Q_q0
            , d.Q_q1, d.Q_q2, d.Q_q3
            , d.SN, d.PDOP, d.HDOP,
            d.VDOP,datacnt);
}

void write_slidefit_data(struct slideRes r, FILE *fp) {
    fprintf(fp, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
            datacnt, r.slideAcc_x, r.slideAcc_y, r.slideAcc_z, r.slideGyo_x, r.slideGyo_y, r.slideGyo_z,
            r.slideMag_x, r.slideMag_y, r.slideMag_z, r.slideJy_yaw, r.slideGps_yaw, r.slideGps_v);
}

void write_raw_gps(struct Data d, FILE *fp) {
    fprintf(fp, "%d,%f,%f,%f\n", datacnt, d.Longitude, d.Lattitude, d.Height);
}

//-----------------------------deal thread end-----------------------------------------

void quaternionToMahonyR(struct Quaternion Q) {
    extern float q0, q1, q2, q3;
    q0 = (float) Q.q1;
    q1 = (float) Q.q2;
    q2 = (float) Q.q3;
    q3 = (float) Q.q4;

    double q0q0 = q0 * q0;
    double q1q1 = q1 * q1;
    double q2q2 = q2 * q2;
    double q3q3 = q3 * q3;
    double q0q1 = q0 * q1;
    double q0q2 = q0 * q2;
    double q0q3 = q0 * q3;
    double q1q2 = q1 * q2;
    double q1q3 = q1 * q3;
    double q2q3 = q2 * q3;
    //定义Cnb
    mahonyR[0] = q0q0 + q1q1 - q2q2 - q3q3;
    mahonyR[1] = 2 * (q1q2 - q0q3);
    mahonyR[2] = 2 * (q1q3 + q0q2);
    mahonyR[3] = 2 * (q1q2 + q0q3);
    mahonyR[4] = q0q0 - q1q1 + q2q2 - q3q3;
    mahonyR[5] = 2 * (q2q3 - q0q1);
    mahonyR[6] = 2 * (q1q3 - q0q2);
    mahonyR[7] = 2 * (q2q3 + q0q1);
    mahonyR[8] = q0q0 - q1q1 - q2q2 + q3q3;
}

struct Quaternion euarToQuaternion(double yaw, double pitch, double roll) {
    struct Quaternion Qresult;
    double cosRoll = cos(roll * 0.5);
    double sinRoll = sin(roll * 0.5);
    double cosPitch = cos(pitch * 0.5);
    double sinPitch = sin(pitch * 0.5);
    double cosYaw = cos(yaw * 0.5);
    double sinYaw = sin(yaw * 0.5);

    Qresult.q1 = (cosRoll * cosPitch * cosYaw + sinRoll * sinPitch * sinYaw);
    Qresult.q2 = (sinRoll * cosPitch * cosYaw - cosRoll * sinPitch * sinYaw);
    Qresult.q3 = (cosRoll * sinPitch * cosYaw + sinRoll * cosPitch * sinYaw);
    Qresult.q4 = (cosRoll * cosPitch * sinYaw - sinRoll * sinRoll * cosYaw);

    return Qresult;
}

int openUART(int fd) {
//----------------------serial port open-----------------------------------------------
    if ((fd = open("/dev/ttyUSB0", O_RDWR | O_NOCTTY)) < 0) {
        printf("open failed\n");
        return -1;
    }
    printf("fd=%d\n", fd);
////////////////setup uart///////////////////////////////
    struct termios newtio, oldtio;
    if (tcgetattr(fd, &oldtio) != 0) {
        printf("get old error\n");
        close(fd);
        return -1;
    }

    bzero(&newtio, sizeof(newtio));
    newtio.c_cflag |= CLOCAL | CREAD;
    newtio.c_cflag &= ~CSIZE;
    newtio.c_cflag |= CS8;
    newtio.c_cflag &= ~PARENB;
    cfsetispeed(&newtio, B9600);
    cfsetospeed(&newtio, B9600);
    newtio.c_cflag &= ~CSTOPB;
    newtio.c_cc[VTIME] = 0;
    newtio.c_cc[VMIN] = 1;
    tcflush(fd, TCIFLUSH);
    if ((tcsetattr(fd, TCSANOW, &newtio)) != 0) {
        printf("set new error\n");
        close(fd);
        return -1;
    }
    return fd;
}

void collectSensorData(int fd) {
    //int fd = *(int *) ft;
    // 50 Time input_data :  51 acc_(xyz) 52 gyo_(xyz) 53 angle_(xyz) 54 m_(xyz) 57 Longitude Lattitude 58 Height Yaw GPSV 59 Q_(q0q1q2q3) 5A (P H V)DOP
    unsigned char buff[2];
    unsigned char ucRxBuffer[250];
    int ucRxCnt = 0;
    int AngCalCnt = 0;
    int AngCalStatus = 0;
    int AngCalN = accCalNum;
    double acc0_x = 0, acc0_y = 0, acc0_z = 0;
    double mag0_x = 0, mag0_y = 0, mag0_z = 0;

    printf("Please wait Eular static calculation..\n");
    //start = clock();
    while (1) {
        read(fd, buff, 1);
        ucRxBuffer[ucRxCnt++] = buff[0];
        //printf("%02x ", buff[i]);
        if (ucRxBuffer[0] != 0x55) {
            ucRxCnt = 0;
            continue;
        }
        if (ucRxCnt < 11) {
//            printf("not enough\n");
            continue;
        } else {
            switch (ucRxBuffer[1]) {
                case 0x50:
                    stcTime.ucYear = ucRxBuffer[2];
                    stcTime.ucMonth = ucRxBuffer[3];
                    stcTime.ucDay = ucRxBuffer[4];
                    stcTime.ucHour = ucRxBuffer[5];
                    stcTime.ucMinute = ucRxBuffer[6];
                    stcTime.ucSecond = ucRxBuffer[7];
                    stcTime.ucMiliSecond = (ucRxBuffer[9] << 8) | ucRxBuffer[8];
                    input_data.Minute = stcTime.ucMinute;
                    input_data.Second = (double) stcTime.ucSecond + (double) stcTime.ucMiliSecond / 1000;
                    break;
                case 0x51:
                    stcAcc.a[0] = ((short) ucRxBuffer[3] << 8) | ucRxBuffer[2];
                    stcAcc.a[1] = ((short) ucRxBuffer[5] << 8) | ucRxBuffer[4];
                    stcAcc.a[2] = ((short) ucRxBuffer[7] << 8) | ucRxBuffer[6];

                    input_data.acc_x = (double) stcAcc.a[0] / 32768 * 16;
                    input_data.acc_y = (double) stcAcc.a[1] / 32768 * 16;
                    input_data.acc_z = (double) stcAcc.a[2] / 32768 * 16;

                    break;
                case 0x52:
                    stcGyro.w[0] = ((short) ucRxBuffer[3] << 8) | ucRxBuffer[2];
                    stcGyro.w[1] = ((short) ucRxBuffer[5] << 8) | ucRxBuffer[4];
                    stcGyro.w[2] = ((short) ucRxBuffer[7] << 8) | ucRxBuffer[6];

                    input_data.gyo_x = (double) stcGyro.w[0] / 32768 * 2000;
                    input_data.gyo_y = (double) stcGyro.w[1] / 32768 * 2000;
                    input_data.gyo_z = (double) stcGyro.w[2] / 32768 * 2000;
                    break;
                case 0x53:
                    stcAngle.Angle[0] = ((short) ucRxBuffer[3] << 8) | ucRxBuffer[2];
                    stcAngle.Angle[1] = ((short) ucRxBuffer[5] << 8) | ucRxBuffer[4];
                    stcAngle.Angle[2] = ((short) ucRxBuffer[7] << 8) | ucRxBuffer[6];

                    input_data.angle_x = (double) stcAngle.Angle[0] / 32768 * 180;
                    input_data.angle_y = (double) stcAngle.Angle[1] / 32768 * 180;
                    input_data.angle_z = (double) stcAngle.Angle[2] / 32768 * 180;
                    break;
                case 0x54:
                    stcMag.h[0] = ((short) ucRxBuffer[3] << 8) | ucRxBuffer[2];
                    stcMag.h[1] = ((short) ucRxBuffer[5] << 8) | ucRxBuffer[4];
                    stcMag.h[2] = ((short) ucRxBuffer[7] << 8) | ucRxBuffer[6];


                    // if (AngCalStatus == 0 && AngCalCnt < AngCalN) {
                        // acc0_x += (double) stcAcc.a[0] / 32768 * 16;
                        // acc0_y += (double) stcAcc.a[1] / 32768 * 16;
                        // acc0_z += (double) stcAcc.a[2] / 32768 * 16;
                        // mag0_x += stcMag.h[0];
                        // mag0_y += stcMag.h[1];
                        // mag0_z += stcMag.h[2];

                        // AngCalCnt++;

                        // printf("initial number is:%d\n", AngCalCnt);
                        // continue;
                    // } else if (AngCalStatus == 0) {
                        // acc0_x /= AngCalN;
                        // acc0_y /= AngCalN;
                        // acc0_z /= AngCalN;
                        // mag0_x /= AngCalN;
                        // mag0_y /= AngCalN;
                        // mag0_z /= AngCalN;

                        // accCalStc_X = acc0_x;
                        // accCalStc_Y = acc0_y;
                        // accCalStc_Z = acc0_z;

                    // } else

					input_data.m_x = stcMag.h[0];
					input_data.m_y = stcMag.h[1];
					input_data.m_z = stcMag.h[2];

                    break;
                case 0x57:
                    stcLonLat.lLon = ((short) ucRxBuffer[5] << 24) | ((short) ucRxBuffer[4] << 16) |
                                     ((short) ucRxBuffer[3] << 8) | ucRxBuffer[2];
                    stcLonLat.lLat = ((short) ucRxBuffer[9] << 24) | ((short) ucRxBuffer[8] << 16) |
                                     ((short) ucRxBuffer[7] << 8) | ucRxBuffer[6];

                    input_data.Longitude = stcLonLat.lLon / 10000000 + (double) (stcLonLat.lLon % 10000000) / 1e5 / 60;
                    input_data.Lattitude = stcLonLat.lLat / 10000000 + (double) (stcLonLat.lLat % 10000000) / 1e5 / 60;

                    break;
                case 0x58:
                    stcGPSV.sGPSHeight = ((short) ucRxBuffer[3] << 8) | ucRxBuffer[2];
                    stcGPSV.sGPSYaw = ((short) ucRxBuffer[5] << 8) | ucRxBuffer[4];
                    stcGPSV.lGPSVelocity = ((short) ucRxBuffer[9] << 24) | ((short) ucRxBuffer[8] << 16) |
                                           ((short) ucRxBuffer[7] << 8) | ucRxBuffer[6];

                    input_data.Height = (double) stcGPSV.sGPSHeight / 10;
                    input_data.Yaw = (double) stcGPSV.sGPSYaw / 100;
                    input_data.GPSV = (double) stcGPSV.lGPSVelocity / 1000;
                    break;
                case 0x59:
                    stcQUAT.sq0 = ((short) ucRxBuffer[3] << 8) | ucRxBuffer[2];
                    stcQUAT.sq1 = ((short) ucRxBuffer[5] << 8) | ucRxBuffer[4];
                    stcQUAT.sq2 = ((short) ucRxBuffer[7] << 8) | ucRxBuffer[6];
                    stcQUAT.sq3 = ((short) ucRxBuffer[9] << 8) | ucRxBuffer[8];

                    input_data.Q_q0 = (double) stcQUAT.sq0 / 32768;
                    input_data.Q_q1 = (double) stcQUAT.sq1 / 32768;
                    input_data.Q_q2 = (double) stcQUAT.sq2 / 32768;
                    input_data.Q_q3 = (double) stcQUAT.sq3 / 32768;
                    break;
                case 0x5A:
                    stcDOP.SN = ((short) ucRxBuffer[3] << 8) | ucRxBuffer[2];
                    stcDOP.sPDOP = ((short) ucRxBuffer[5] << 8) | ucRxBuffer[4];
                    stcDOP.sHDOP = ((short) ucRxBuffer[7] << 8) | ucRxBuffer[6];
                    stcDOP.sVDOP = ((short) ucRxBuffer[9] << 8) | ucRxBuffer[8];

                    input_data.SN = stcDOP.SN;
                    input_data.PDOP = (double) stcDOP.sPDOP / 100;
                    input_data.HDOP = (double) stcDOP.sHDOP / 100;
                    input_data.VDOP = (double) stcDOP.sVDOP / 100;
                    EnQueue(&queue, input_data);
                    break;
            }


        }
        ucRxCnt = 0;
        if ((queue.rear - queue.front) == 1) {
            printf("The Acc static drift is: \n driftX = %.5f, driftY = %.5f, driftZ = %.5f\r\n", (float) accCalStc_X,
                   (float) accCalStc_Y, (float) accCalStc_Z);
            printf("The first Eular angles is: \n pitch0 = %.5f, roll0 = %.5f, yaw0 = %.5f\r\n", (float) pitch0,
                   (float) roll0, (float) yaw0);
            printf("Time:20%d-%d-%d %d:%d:%.3f\r\n", stcTime.ucYear, stcTime.ucMonth, stcTime.ucDay, stcTime.ucHour,
                   stcTime.ucMinute, (float) stcTime.ucSecond + (float) stcTime.ucMiliSecond / 1000);
            printf("Acc:%.3f %.3f %.3f\r\n", (float) stcAcc.a[0] / 32768 * 16 - accCalStc_X,
                   (float) stcAcc.a[1] / 32768 * 16 - accCalStc_Y, (float) stcAcc.a[2] / 32768 * 16 - accCalStc_Z);
            printf("Gyro:%.3f %.3f %.3f\r\n", (float) stcGyro.w[0] / 32768 * 2000, (float) stcGyro.w[1] / 32768 * 2000,
                   (float) stcGyro.w[2] / 32768 * 2000);
            printf("Angle:%.3f %.3f %.3f\r\n", (float) stcAngle.Angle[0] / 32768 * 180,
                   (float) stcAngle.Angle[1] / 32768 * 180, (float) stcAngle.Angle[2] / 32768 * 180);
            printf("Mag:%d %d %d\r\n", stcMag.h[0], stcMag.h[1], stcMag.h[2]);
            printf("Longitude:%.8f Lattitude:%.8f\r\n",
                   stcLonLat.lLon / 10000000 + (double) (stcLonLat.lLon % 10000000) / 1e5 / 60,
                   stcLonLat.lLat / 10000000 + (double) (stcLonLat.lLat % 10000000) / 1e5 / 60);
            printf("GPSHeight:%.1fm GPSYaw:%.3fDeg GPSV:%.3fkm/h\r\n", (float) stcGPSV.sGPSHeight / 10,
                   (float) stcGPSV.sGPSYaw / 100, (float) stcGPSV.lGPSVelocity / 1000);
//			printf("Quaternion:%.3f %.3f %.3f %.3f\r\n",(float)stcQUAT.sq0/32768,(float)stcQUAT.sq1/32768,(float)stcQUAT.sq2/32768,(float)stcQUAT.sq3/32768);
            printf("Quaternion:%.3f %.3f %.3f %.3f\r\n", (float) Q0.q1, (float) Q0.q2, (float) Q0.q3, (float) Q0.q4);
            printf("GPS number:%d GPSDOP:PDOP%.7f HDOP%.7f VDOP%.7f\r\n", stcDOP.SN, (float) stcDOP.sPDOP / 100,
                   (float) stcDOP.sHDOP / 100, (float) stcDOP.sVDOP / 100);

        }
    }
}

//------------------------------main thread--------------------------------------------
int main(int argc, char *argv[]) {
    smoothRes = (double ***) malloc(sizeof(double **) * 4);
    for (int i = 0; i < 4; i++) {
        smoothRes[i] = (double **) malloc(sizeof(double *) * 3);
        for (int j = 0; j < 3; j++)
            smoothRes[i][j] = (double *) malloc(sizeof(double) * 12);
    }


    int fd = 0;
    mahonyR = (double *) malloc(sizeof(double) * 9);
    Fn = (double *) malloc(sizeof(double) * 3);
    //加速度常值飘逸
    //accStaticDriftCal();
    //printf("1");
    pthread_t thread1, thread2;
    pthread_mutex_t mutex;
    fd = openUART(fd);
    if (fd == -1)
        return -1;
    pthread_create(&thread1, NULL, start_deal, &queue);
//    pthread_create(&thread2, NULL, collectSensorData, &fd);

    collectSensorData(fd);
    close(fd);
    return 0;
}

