%SinIn[i] = A0 * (sin(2*PI*i/25)+sin(2*PI* i * 0.4 ));
clc
clear
N=128;A0=255;PI=3.1415926;
i=1:N;
x=A0 * (sin(2*PI*i/25)+sin(2*PI* i * 0.4 ));
y=fft(x,N);%���и���Ҷ�任��÷�����
z=abs(y);%����õķ����׼Ӿ���ֵ   
k=angle(y);%����λ��
subplot(2,3,1);plot(i,x);title('MATLAB�е����벨��');
subplot(2,3,2);plot(i,z);title('MATLAB�ķ�����');
subplot(2,3,3);plot(i,k);title('MATLAB����λ��');
var1=load('input.txt');%��ȡC������N����������
var2=load('output.txt');%��ȡC������N��FFT��������
var3=load('outputt.txt');%��ȡC������N����λ����
subplot(2,3,4);plot(i,var1);title('C�е����벨��');
subplot(2,3,5);plot(i,abs(var2));title('C�ķ�����');
subplot(2,3,6);plot(i,abs(var3));title('C����λ��');



