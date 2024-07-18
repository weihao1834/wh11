#include <math.h>
#include <stdlib.h>
#include <conio.h>
#include <stdio.h>
#define PI 3.1415926535
#define N 128
#define M 7
#define A0 255.0 //�������벨�εķ�ֵ

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
	//���벨��
	for (i = 0; i < N; i++)
	{
		SinIn[i] = A0 * (sin(2 * PI * i / 25) + sin(2 * PI * i * 0.4));
		dataR[i] = SinIn[i]; //���벨�ε�ʵ������
		dataI[i] = 0;		 //���벨�ε���������
		dataA[i] = 0;		 //������εķ�����Ϊ0
	}
	if ((fp = fopen("input.txt", "wb")) == NULL)
	{
		printf("�޷����ļ�");
		exit(0); //��ֹ����
	}
	for (i = 0; i < N; i++)
		fprintf(fp, "%f ", dataR[i]);
	//fwrite(&dataR,sizeof(dataR),1,fp);
	fclose(fp);
	//�������е���
	/*
J=N/2;
for(I = 1; I < N - 2; I++)
{
	//����Ԫ�ػ���
	if(I<J)
	{
		temp = dataR[I];
		dataR[I] = dataR[J];
		dataR[J] = temp;
	}
	//��ȡ��һ��������
	K=N/2;//K��N��Ȩֵ�����������λȨֵ��ʼ�жϡ�
	while(1)
	{
		//������ѭ���жϸ�λ�Ƿ�Ϊ1
 		if(J< K)	
		{
			J=J+K;//�жϽ��Ϊ��λΪ0�����ϸ�λȨֵ��0���1��
			break;//ѭ������
 		}
		else 
		{
			J=J-K;//�жϽ��Ϊ��λΪ1����ȥ��λȨֵ��1���0��
			K=K/2;//�õ���һλ��Ȩֵ��������һλ0��1���жϡ�
		 }
	}
}*/

	//����
	for (I = 0; I < N; I++) //���ݹ����ģ���Ҫ������Ԫ��ִ����䵹��
	{
		/*��ȡ�±�I�ķ���J����ֵ*/
		J = 0;
		for (k = 0; k < (M / 2 + 0.5); k++) //k��ʾ����
		{
			//*�������*/
			m = 1;							 //m�����λΪ1�Ķ�������
			n = (unsigned int)pow(2, M - 1); //n�ǵ�MλΪ1�Ķ�������
			m <<= k;						 //��m����kλ
			n >>= k;						 //��n����kλ
			F0 = I & n;						 //I��n��λ����ȡ��ǰ�벿�ֵ�kλ
			F1 = I & m;						 //I��m��λ����ȡ��F0��Ӧ�ĺ�벿�ֵĵ�λ
			if (F0)
				J = J | m; //J��m��λ��ʹF0��Ӧ��λΪ1
			if (F1)
				J = J | n; //J��n��λ��ʹF1��Ӧ��λΪ1
		}
		//printf("IΪ��");printBin(I,M);printf(" %d\n",I);
		//printf("JΪ��");printBin(J,M);printf(" %d\n\n",J);

		if (I < J)
		{
			temp = dataR[I];
			dataR[I] = dataR[J];
			dataR[J] = temp;
		}
	}

	//����FFT
	/*	Tr = Xr(J+B)cos(2.0*PI*p/N) + Xi(J+B)sin(2.0*PI*p/N)
	    Ti = Xi(J+B)cos(2.0*PI*p/N) - Xr(J+B)sin(2.0*PI*p/N)
		Ar[J] = Xr(J) + Tr
	    Ai[J] = Xi(J) + Ti
	    Ar[J+B] = Xr(J) - Tr
	    Ai[J+B] = Xi(J) - Ti
	   (���� XrΪ��һ����Ar, XiΪ��һ����Ai)*/
	for (L = 1; L <= M; L++) //FFT���μ���L��1--M
	{
		/*��L��������:
		Ȼ����ڵ�L�������������������ᵽ�����������������Ŀ���ڼ��B:�ж����ֵ��������Ҫ��Ҫѭ�����ٴ�;
		����ѭ���Ĳ�ͬ����תָ��PҲ��ͬ��P������Ϊk=2^(M-L)*/
		//�ȼ���һ�¼�� B = 2^(L-1);
		B = 1;
		B = (int)pow(2, L - 1);
		for (j = 0; j <= B - 1; j++)
		//j = 0,1,2,...,2^(L-1) - 1
		{ /*ͬ�ֵ�������*/
			//�ȼ�������k=2^(M-L)
			k = 1;
			k = (int)pow(2, M - L);
			//������תָ��p������Ϊkʱ����P=j*k
			p = 1;
			p = j * k;
			/*�������������������ǿ���֪��ͬ�ֵ�������Ĵ����պ�Ϊ����k=2^(M-L)��
			ͬ�ֵ��ε����������Ϊ��������Ĵ���*/
			for (i = 0; i <= k - 1; i++)
			{
				//�����±궨Ϊr
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

	//�������Ƶ����,dataA[i] = sqrt(dataR[i]*dataR[i]+dataI[i]*dataI[i])
	for (i = 0; i < N; i++)
	{
		dataA[i] = sqrt(dataR[i] * dataR[i] + dataI[i] * dataI[i]);
	}
	if ((fp = fopen("output.txt", "wb")) == NULL)
	{
		printf("�޷����ļ�");
		exit(0); //��ֹ����
	}
	for (i = 0; i < N; i++)
		fprintf(fp, "%f\n", dataA[i]);
	fclose(fp);
	//�����λ�ף���[i] = arctan(dataI[i]  /  dataR[i])
	for (i = 0; i < N; i++)
	{
		dataP[i] = atan(dataI[i] / dataR[i]);
	}
	if ((fp = fopen("outputt.txt", "wb")) == NULL)
	{
		printf("�޷����ļ�");
		exit(0); //��ֹ����
	}
	for (i = 0; i < N; i++)
		fprintf(fp, "%f\n", dataP[i]);
	//printf("%f\n",dataA[i]);}
	//printf("%d\n",M/2);}
	//fwrite(&dataR,sizeof(dataR),1,fp);
	fclose(fp);
	//Ƶ�׵�ƽ����Ϊ�����ף���Ϊ��P[i] = dataR[i]*dataR[i]+dataI[i]*dataI[i]
}
