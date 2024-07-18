%SinIn[i] = A0 * (sin(2*PI*i/25)+sin(2*PI* i * 0.4 ));
clc
clear
N=128;A0=255;PI=3.1415926;
i=1:N;
x=A0 * (sin(2*PI*i/25)+sin(2*PI* i * 0.4 ));
y=fft(x,N);%进行傅里叶变换求得幅度谱
z=abs(y);%对求得的幅度谱加绝对值   
k=angle(y);%求相位谱
subplot(2,3,1);plot(i,x);title('MATLAB中的输入波形');
subplot(2,3,2);plot(i,z);title('MATLAB的幅度谱');
subplot(2,3,3);plot(i,k);title('MATLAB的相位谱');
var1=load('input.txt');%读取C产生的N点序列数据
var2=load('output.txt');%读取C产生的N点FFT序列数据
var3=load('outputt.txt');%读取C产生的N点相位数据
subplot(2,3,4);plot(i,var1);title('C中的输入波形');
subplot(2,3,5);plot(i,abs(var2));title('C的幅度谱');
subplot(2,3,6);plot(i,abs(var3));title('C的相位谱');



