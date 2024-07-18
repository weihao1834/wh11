#include <math.h>
#include <stdlib.h>
#include <conio.h>
#include <stdio.h>
#define PI 3.1415926535
#define N 128
#define M 7
#define A0 255.0 //定义输入波形的幅值

void printBin(unsigned i, int bits)
{
	unsigned shiftcount = bits;
	while (bits-- != 0)
		i &(1 << bits) ? printf("1") : printf("0");
}
void main(void)
{
	int i, j, k, r;
	int p, L, B;
	unsigned int I, J, K, F0, F1, m, n;
	float SinIn[N];
	float dataR[N], dataI[N];
	float dataA[N], dataP[N];
	float Tr, Ti, temp;
	FILE *fp;
	//输入波形
	for (i = 0; i < N; i++)
	{
		SinIn[i] = A0 * (sin(2 * PI * i / 25) + sin(2 * PI * i * 0.4));
		dataR[i] = SinIn[i]; //输入波形的实数部分
		dataI[i] = 0;		 //输入波形的虚数部分
		dataA[i] = 0;		 //输出波形的幅度谱为0
	}
	if ((fp = fopen("input.txt", "wb")) == NULL)
	{
		printf("无法打开文件");
		exit(0); //终止程序
	}
	for (i = 0; i < N; i++)
		fprintf(fp, "%f ", dataR[i]);
	//fwrite(&dataR,sizeof(dataR),1,fp);
	fclose(fp);
	//输入序列倒序
	/*
J=N/2;
for(I = 1; I < N - 2; I++)
{
	//数组元素互换
	if(I<J)
	{
		temp = dataR[I];
		dataR[I] = dataR[J];
		dataR[J] = temp;
	}
	//求取下一个反序数
	K=N/2;//K是N的权值变量，从最高位权值开始判断。
	while(1)
	{
		//从左到右循环判断该位是否为1
 		if(J< K)	
		{
			J=J+K;//判断结果为该位为0，加上该位权值，0变成1。
			break;//循环结束
 		}
		else 
		{
			J=J-K;//判断结果为该位为1，减去该位权值，1变成0。
			K=K/2;//得到下一位的权值，继续下一位0或1的判断。
		 }
	}
}*/

	//倒序
	for (I = 0; I < N; I++) //根据规律四，需要对数组元素执行码间倒序
	{
		/*获取下标I的反序J的数值*/
		J = 0;
		for (k = 0; k < (M / 2 + 0.5); k++) //k表示操作
		{
			//*反序操作*/
			m = 1;							 //m是最低位为1的二进制数
			n = (unsigned int)pow(2, M - 1); //n是第M位为1的二进制数
			m <<= k;						 //对m左移k位
			n >>= k;						 //对n右移k位
			F0 = I & n;						 //I与n按位与提取出前半部分第k位
			F1 = I & m;						 //I与m按位与提取出F0对应的后半部分的低位
			if (F0)
				J = J | m; //J与m按位或使F0对应低位为1
			if (F1)
				J = J | n; //J与n按位或使F1对应高位为1
		}
		//printf("I为：");printBin(I,M);printf(" %d\n",I);
		//printf("J为：");printBin(J,M);printf(" %d\n\n",J);

		if (I < J)
		{
			temp = dataR[I];
			dataR[I] = dataR[J];
			dataR[J] = temp;
		}
	}

	//进行FFT
	/*	Tr = Xr(J+B)cos(2.0*PI*p/N) + Xi(J+B)sin(2.0*PI*p/N)
	    Ti = Xi(J+B)cos(2.0*PI*p/N) - Xr(J+B)sin(2.0*PI*p/N)
		Ar[J] = Xr(J) + Tr
	    Ai[J] = Xi(J) + Ti
	    Ar[J+B] = Xr(J) - Tr
	    Ai[J+B] = Xi(J) - Ti
	   (其中 Xr为上一级的Ar, Xi为上一级的Ai)*/
	for (L = 1; L <= M; L++) //FFT蝶形级数L从1--M
	{
		/*第L级的运算:
		然后对于第L级，我们在问题五种提到，蝶形运算的种类数目等于间隔B:有多少种蝶形运算就要需要循环多少次;
		随着循环的不同，旋转指数P也不同，P的增量为k=2^(M-L)*/
		//先计算一下间隔 B = 2^(L-1);
		B = 1;
		B = (int)pow(2, L - 1);
		for (j = 0; j <= B - 1; j++)
		//j = 0,1,2,...,2^(L-1) - 1
		{ /*同种蝶形运算*/
			//先计算增量k=2^(M-L)
			k = 1;
			k = (int)pow(2, M - L);
			//计算旋转指数p，增量为k时，则P=j*k
			p = 1;
			p = j * k;
			/*接下来，由问题六我们可以知道同种蝶形运算的次数刚好为增量k=2^(M-L)；
			同种蝶形的运算次数即为蝶形运算的次数*/
			for (i = 0; i <= k - 1; i++)
			{
				//数组下标定为r
				r = 1;
				r = j + 2 * B * i;
				Tr = dataR[r + B] * cos(2.0 * PI * p / N) + dataI[r + B] * sin(2.0 * PI * p / N);
				Ti = dataI[r + B] * cos(2.0 * PI * p / N) - dataR[r + B] * sin(2.0 * PI * p / N);
				dataR[r + B] = dataR[r] - Tr;
				dataI[r + B] = dataI[r] - Ti;
				dataR[r] = dataR[r] + Tr;
				dataI[r] = dataI[r] + Ti;
			}
		}
	}

	//求出幅度频率谱,dataA[i] = sqrt(dataR[i]*dataR[i]+dataI[i]*dataI[i])
	for (i = 0; i < N; i++)
	{
		dataA[i] = sqrt(dataR[i] * dataR[i] + dataI[i] * dataI[i]);
	}
	if ((fp = fopen("output.txt", "wb")) == NULL)
	{
		printf("无法打开文件");
		exit(0); //终止程序
	}
	for (i = 0; i < N; i++)
		fprintf(fp, "%f\n", dataA[i]);
	fclose(fp);
	//求出相位谱，Φ[i] = arctan(dataI[i]  /  dataR[i])
	for (i = 0; i < N; i++)
	{
		dataP[i] = atan(dataI[i] / dataR[i]);
	}
	if ((fp = fopen("outputt.txt", "wb")) == NULL)
	{
		printf("无法打开文件");
		exit(0); //终止程序
	}
	for (i = 0; i < N; i++)
		fprintf(fp, "%f\n", dataP[i]);
	//printf("%f\n",dataA[i]);}
	//printf("%d\n",M/2);}
	//fwrite(&dataR,sizeof(dataR),1,fp);
	fclose(fp);
	//频谱的平方称为功率谱，记为：P[i] = dataR[i]*dataR[i]+dataI[i]*dataI[i]
}
