#include <iostream.h>
#include <stdio.h>
#include <math.h>
double T;//定义时间t的最大范围
#define N 2000
#define e 0.0001
#define r 1.4
static double U[N+1];
static double U1[N+1],U2[N+1],U3[N+1];
static double u[N+1],d[N+1],p[N+1];
double dx=1.0/N;
double ffd[N+1];
double fzd[N+1];
double *FZ(double fz[N+1])//定义a>0时的WENO差分函数
{
	double m=0.000001;
    double IS1z[N+1],IS2z[N+1],IS3z[N+1];
	int i;
	for(i=0;i<N+1;i++)//生成IS1[N+1]数组
	{
		if(i>1&&i<N+1)
		{
			IS1z[i]=0.25*(fz[i-2]-4*fz[i-1]+3*fz[i])*(fz[i-2]-4*fz[i-1]+3*fz[i])+13.0/12*(fz[i-2]-2*fz[i-1]+fz[i])*(fz[i-2]-2*fz[i-1]+fz[i]);
		}
		if(i==0)
		{
			IS1z[i]=0;
		}
		if(i==1)
		{
			IS1z[i]=10.0/3*(fz[i]-fz[i-1])*(fz[i]-fz[i-1]);
		}
	}
	for(i=0;i<N+1;i++)//生成IS2[N+1]数组
	{
		if(i>0&&i<N)
		{
			IS2z[i]=1.0/4*(fz[i-1]-fz[i+1])*(fz[i-1]-fz[i+1])+13.0/12*(fz[i-1]-2*fz[i]+fz[i+1])*(fz[i-1]-2*fz[i]+fz[i+1]);
		}
		if(i==0)
		{
			IS2z[i]=0;
		}
		if(i==N)
		{
			IS2z[i]=4.0/3*(fz[i-1]-fz[i])*(fz[i-1]-fz[i]);
		}
	}
	for(i=0;i<N+1;i++)//生成IS3[N+1]数组
	{
		if(i>=0&&i<N-1)
		{
			IS3z[i]=1.0/4*(3*fz[i]-4*fz[i+1]+fz[i+2])*(3*fz[i]-4*fz[i+1]+fz[i+2])+13.0/12*(fz[i]-2*fz[i+1]+fz[i+2])*(fz[i]-2*fz[i+1]+fz[i+2]);
		}
		if(i==N-1)
		{
			IS3z[i]=10.0/3*(fz[i]-fz[i+1])*(fz[i]-fz[i+1]);
		}
		if(i==N)
		{
			IS3z[i]=0;
		}
	}
	double Alph1z[N+1],Alph2z[N+1],Alph3z[N+1];//定义三个Alph数组
	for(i=0;i<N+1;i++)//生成3个Alph数组
	{
		Alph1z[i]=1.0/10/(m+IS1z[i])/(m+IS1z[i]);
        Alph2z[i]=6.0/10/(m+IS2z[i])/(m+IS2z[i]);
		Alph3z[i]=3.0/10/(m+IS3z[i])/(m+IS3z[i]);
	}
	double w1z[N+1],w2z[N+1],w3z[N+1];
    for(i=0;i<N+1;i++)
	{
		if(i>1&&i<N+1)
		{
		w1z[i]=Alph1z[i]*1.0/(Alph1z[i]+Alph2z[i]+Alph3z[i]);
		}
		else
		{
		w1z[i]=0;
		}
	}
	for(i=0;i<N+1;i++)
	{
		if(i>0&&i<N)
		{
        w2z[i]=Alph2z[i]*1.0/(Alph1z[i]+Alph2z[i]+Alph3z[i]);
		}
		else
		{
		w2z[i]=0;
		}
	}
	for(i=0;i<N+1;i++)
	{
		if(i>=0&&i<N-1)
		{
		w3z[i]=Alph3z[i]*1.0/(Alph1z[i]+Alph2z[i]+Alph3z[i]);
		}
		else
		{
		w3z[i]=0;
		}
	}
	double fz1[N+1],fz2[N+1],fz3[N+1];
	for(i=0;i<N+1;i++)//通过这个循环产生fz1[N+1]数组的元素
	{
		if(i>1&&i<N+1)
		{
			fz1[i]=1.0/3*fz[i-2]-7.0/6*fz[i-1]+11.0/6*fz[i];
		}
		else
		{
            fz1[i]=0;
		}
	}
	for(i=0;i<N+1;i++)//通过这个循环产生fz2[N+1]数组的元素
	{
		if(i>0&&i<N)
		{
			fz2[i]=-1.0/6*fz[i-1]+5.0/6*fz[i]+1.0/3*fz[i+1];
		}
        else
		{
			fz2[i]=0;
		}
	}
	for(i=0;i<N+1;i++)//通过这个循环产生fz3[N+1]数组的元素
	{
		if(i>=0&&i<N-1)
		{
			fz3[i]=1.0/3*fz[i]+5.0/6*fz[i+1]-1.0/6*fz[i+2];
		}
        else
		{
			fz3[i]=0;
		}
	}//至此，3个数组生成完毕
	double fzweno[N+1];
	for(i=0;i<N+1;i++)//生成fzweno[N+1]数组的元素
	{
		fzweno[i]=w1z[i]*fz1[i]+w2z[i]*fz2[i]+w3z[i]*fz3[i];
	}
	for(i=0;i<N+1;i++)
	{
		if(i<=2||i>=N-2)
		{
			fzd[i]=0;
		}
		if(i>2&&i<N-2)
		{
			fzd[i]=(fzweno[i]-fzweno[i-1])*1.0/dx;
		}
	}
	return fzd;
}
double *FF(double ff[N+1])//定义a<0时的WENO差分函数
{
	double n=0.000001;
    double IS1f[N+1],IS2f[N+1],IS3f[N+1];
	int i;
	for(i=0;i<N+1;i++)//生成IS1[N+1]数组
	{
		if(i>=0&&i<N-1)
		{
			IS1f[i]=1.0/4*(ff[i+2]-4*ff[i+1]+3*ff[i])*(ff[i+2]-4*ff[i+1]+3*ff[i])+13.0/12*(ff[i+2]-2*ff[i+1]+ff[i])*(ff[i+2]-2*ff[i+1]+ff[i]);
		}
		if(i==N)
		{
			IS1f[i]=0;
		}
		if(i==N-1)
		{
			IS1f[i]=10.0/3*(ff[i]-ff[i-1])*(ff[i]-ff[i-1]);
		}
	}
	for(i=0;i<N+1;i++)//生成IS2[N+1]数组
	{
		if(i>0&&i<N)
		{
			IS2f[i]=1.0/4*(ff[i+1]-ff[i-1])*(ff[i+1]-ff[i-1])+13.0/12*(ff[i+1]-2*ff[i]+ff[i-1])*(ff[i+1]-2*ff[i]+ff[i-1]);
		}
		if(i==0)
		{
			IS2f[i]=0;
		}
		if(i==N)
		{
			IS2f[i]=4.0/3*(ff[i-1]-ff[i])*(ff[i-1]-ff[i]);
		}
	}
	for(i=0;i<N+1;i++)//生成IS3[N+1]数组
	{
		if(i>1&&i<N+1)
		{
			IS3f[i]=1.0/4*(3*ff[i]-4*ff[i-1]+ff[i-2])*(3*ff[i]-4*ff[i-1]+ff[i-2])+13.0/12*(ff[i]-2*ff[i-1]+ff[i-2])*(ff[i]-2*ff[i-1]+ff[i-2]);
		}
		if(i==1)
		{
			IS3f[i]=10.0/3*(ff[i]-ff[i-1])*(ff[i]-ff[i-1]);
		}
		if(i==0)
		{
			IS3f[i]=0;
		}
	}
	double Alph1f[N+1],Alph2f[N+1],Alph3f[N+1];//定义三个Alph数组
	for(i=0;i<N+1;i++)//生成3个Alph数组
	{
		Alph1f[i]=1.0/10/(n+IS1f[i])/(n+IS1f[i]);
        Alph2f[i]=6.0/10/(n+IS2f[i])/(n+IS2f[i]);
		Alph3f[i]=3.0/10/(n+IS3f[i])/(n+IS3f[i]);
	}
	double w1f[N+1],w2f[N+1],w3f[N+1];
    for(i=0;i<N+1;i++)
	{
		if(i>=0&&i<N-1)
		{
		w1f[i]=Alph1f[i]*1.0/(Alph1f[i]+Alph2f[i]+Alph3f[i]);
		}
		else
		{
		w1f[i]=0;//边界处采用简易的降阶处理方法
		}
	}
    for(i=0;i<N+1;i++)
	{
		if(i>0&&i<N)
		{
        w2f[i]=Alph2f[i]*1.0/(Alph1f[i]+Alph2f[i]+Alph3f[i]);
		}
		else
		{
		w2f[i]=0;
		}
	}
	for(i=0;i<N+1;i++)
	{
		if(i>1&&i<N+1)
		{
		w3f[i]=Alph3f[i]*1.0/(Alph1f[i]+Alph2f[i]+Alph3f[i]);
		}
		else
		{
		w3f[i]=0;
		}
	}
	double ff1[N+1],ff2[N+1],ff3[N+1];
	for(i=0;i<N+1;i++)//通过这个循环产生ff1[N+1]数组的元素
	{
		if(i>=0&&i<N-1)
		{
			ff1[i]=1.0/3*ff[i+2]-7.0/6*ff[i+1]+11.0/6*ff[i];
		}
		else
		{
			ff1[i]=0;
		}
	}
	for(i=0;i<N+1;i++)//通过这个循环产生ff2[N+1]数组的元素
	{
		if(i>0&&i<N)
		{
			ff2[i]=-1.0/6*ff[i+1]+5.0/6*ff[i]+1.0/3*ff[i-1];
		}
		else
		{
			ff2[i]=0;
		}
	}
	for(i=0;i<N+1;i++)//通过这个循环产生fz3[N+1]数组的元素
	{
		if(i>1&&i<N+1)
		{
			ff3[i]=1.0/3*ff[i]+5.0/6*ff[i-1]-1.0/6*ff[i-2];
		}
		else
		{
			ff3[i]=0;
		}
	}//至此，3个数组生成完毕
	double ffweno[N+1];
	for(i=0;i<N+1;i++)//生成fzweno[N+1]数组的元素
	{
		ffweno[i]=w1f[i]*ff1[i]+w2f[i]*ff2[i]+w3f[i]*ff3[i];
	}
    double ffd[N+1];
	for(i=0;i<N+1;i++)
	{
		if(i<=2||i>=N-2)
		{
			ffd[i]=0;
		}
		if(i>2&&i<N-2)
		{
			ffd[i]=(ffweno[i+1]-ffweno[i])*1.0/dx;//修改
		}
	}
	return ffd;
}
void main()
{
	double t=0;
	double dt;
	cout<<"请输入求解的时刻："<<endl;
	cin>>T;
	cout<<"您输入的求解的时刻是："<<T<<"时刻"<<endl;
    cout<<"请输入求解的时间步长："<<endl;
	cin>>dt;
	cout<<"您输入的时间步长是："<<dt<<endl;
    int j=0;
	for(j=0;j<N+1;j++)//先通过一个循环结构产生N+1个点处的密度初值
	{
	   if(j<N/10)
	   {
		   d[j]=3.857;
	   }
	   if(j>N/10)
	   {
		   d[j]=1+0.3*sin(40*j*1.0/N);
	   }
	}
	for(j=0;j<N+1;j++)//先通过一个循环结构产生N+1个点处的压力初值
	{
	   if(j<N/10)
	   {
		   p[j]=10.333;
	   }
	   if(j>=N/10)
	   {
		   p[j]=1.0;
	   }
	}
	for(j=0;j<N+1;j++)//先通过一个循环结构产生N+1个点处的速度初值
	{
	   if(j<N/10)
	   {
		   u[j]=2.629;
	   }
	   if(j>=N/10)
	   {
		   u[j]=0;
	   }
	}
		static double u1a[N+1],u1b[N+1],u1c[N+1];//这些分别对应R-K方法中的中间量u1
		static double Lz1[N+1],Lz2[N+1],Lz3[N+1],Lf1[N+1],Lf2[N+1],Lf3[N+1];//定义3阶R-K方法参数
		static double L1[N+1],L2[N+1],L3[N+1];
	while(t<=T)
	{	
		t=t+dt;
	    for(j=0;j<N+1;j++)//通过一个循环结构产生n时刻的U值,方便R-K方法中迭代
		{
            U1[j]=d[j];
			U2[j]=d[j]*u[j];
			U3[j]=p[j]*1.0/(r-1)+0.5*d[j]*u[j]*u[j];
		}  
		static double c[N+1];//定义声速
	    for(j=0;j<N+1;j++)//先通过一个循环结构产生N+1个点处的声速
		{
            c[j]=sqrt(r*p[j]*1.0/d[j]);       
		}
		static double a1[N+1],a2[N+1],a3[N+1];//定义3个特征值
	    for(j=0;j<N+1;j++)//先通过一个循环结构产生N+1个特征值
		{
            a1[j]=u[j];
			a2[j]=u[j]-c[j];
			a3[j]=u[j]+c[j];
		}  
		static double a1z[N+1],a1f[N+1],a2z[N+1],a2f[N+1],a3z[N+1],a3f[N+1];//定义分裂变量
		for(j=0;j<N+1;j++)//为分裂变量赋值
		{
            a1z[j]=(a1[j]+sqrt(a1[j]*a1[j]+e*e))*1.0/2;
            a1f[j]=(a1[j]-sqrt(a1[j]*a1[j]+e*e))*1.0/2;
			a2z[j]=(a2[j]+sqrt(a2[j]*a2[j]+e*e))*1.0/2;
			a2f[j]=(a2[j]-sqrt(a2[j]*a2[j]+e*e))*1.0/2;
			a3z[j]=(a3[j]+sqrt(a3[j]*a3[j]+e*e))*1.0/2;
			a3f[j]=(a3[j]-sqrt(a3[j]*a3[j]+e*e))*1.0/2;
		}
		static double fz1a[N+1],fz2a[N+1],fz3a[N+1];
		static double ff1a[N+1],ff2a[N+1],ff3a[N+1];
		static double wza[N+1],wfa[N+1];
		for(j=0;j<N+1;j++)//为每一个fz及ff赋值
		{
            fz1a[j]=d[j]*1.0/2/r*(2*(r-1)*a1z[j]+a2z[j]+a3z[j]);
            ff1a[j]=d[j]*1.0/2/r*(2*(r-1)*a1f[j]+a2f[j]+a3f[j]);
            fz2a[j]=d[j]*1.0/2/r*(2*(r-1)*a1z[j]*u[j]+a2z[j]*(u[j]-c[j])+a3z[j]*(u[j]+c[j]));
            ff2a[j]=d[j]*1.0/2/r*(2*(r-1)*a1f[j]*u[j]+a2f[j]*(u[j]-c[j])+a3f[j]*(u[j]+c[j]));
            wza[j]=(3-r)*(a2z[j]+a3z[j])*c[j]*c[j]*1.0/2/(r-1);
            wfa[j]=(3-r)*(a2f[j]+a3f[j])*c[j]*c[j]*1.0/2/(r-1);
            fz3a[j]=d[j]*1.0/2/r*((r-1)*a1z[j]*u[j]*u[j]+a2z[j]*(u[j]-c[j])*(u[j]-c[j])*1.0/2+a3z[j]*(u[j]+c[j])*(u[j]+c[j])*1.0/2+wza[j]);
            ff3a[j]=d[j]*1.0/2/r*((r-1)*a1f[j]*u[j]*u[j]+a2f[j]*(u[j]-c[j])*(u[j]-c[j])*1.0/2+a3f[j]*(u[j]+c[j])*(u[j]+c[j])*1.0/2+wfa[j]);
		}
	    int k=0;
		double *fzd=FZ(fz1a);//调用WENO差分函数
		for(k=0;k<N+1;k++)
		{
			Lz1[k]=0-fzd[k];          
		}
		double *ffd=FF(ff1a);
		for(k=0;k<N+1;k++)
		{
			Lf1[k]=0-ffd[k];          
		}
		for(k=0;k<N+1;k++)
		{
			L1[k]=Lz1[k]+Lf1[k];
		}
	    fzd=FZ(fz2a);
		for(k=0;k<N+1;k++)
		{
			Lz2[k]=0-fzd[k];          
		}
		ffd=FF(ff2a);
		for(k=0;k<N+1;k++)
		{
			Lf2[k]=0-ffd[k];          
		}
		for(k=0;k<N+1;k++)
		{
			L2[k]=Lz2[k]+Lf2[k];
		}
		fzd=FZ(fz3a);
		for(k=0;k<N+1;k++)
		{
			Lz3[k]=0-fzd[k];          
		}
		ffd=FF(ff3a);
		for(k=0;k<N+1;k++)
		{
			Lf3[k]=0-ffd[k];          
		}
		for(k=0;k<N+1;k++)
		{
			L3[k]=Lz3[k]+Lf3[k];
		}
		for(k=0;k<N+1;k++)//第一个循环结构，为了求出数组u1[21]
		{
			u1a[k]=U1[k]+dt*L1[k];
			u1b[k]=U2[k]+dt*L2[k];
			u1c[k]=U3[k]+dt*L3[k];
		}//3阶R-K方法第一步结束
	    //求出u1a[k]，u1b[k]，u1c[k]对应的密度，速度，压力数组；再重复求特征值等步骤，按照PPT上的方式，求出u2a[k]，u2b[k]，u2c[k]
	        static double u1[N+1],d1[N+1],p1[N+1];//分别对应R-K方法第一步产生的新的速度、密度、压力
		    for(j=0;j<N+1;j++)//根据u1a[k],u1b[k],u1c[k]产生新的u1[k],d1[k],p1[k];
			{
				d1[j]=u1a[j];
				u1[j]=u1b[j]*1.0/u1a[j];
				p1[j]=(u1c[j]-d1[j]*u1[j]*u1[j]*1.0/2)*(r-1);
			}
        	static double c1[N+1];
			for(j=0;j<N+1;j++)//重新计算声速
			{
                c1[j]=sqrt(r*p1[j]*1.0/d1[j]);//EEEEEEEEEEEEEEE
			}
		    static double b1[N+1],b2[N+1],b3[N+1];//定义3个特征值
	        for(j=0;j<N+1;j++)//先通过一个循环结构产生N+1个特征值
			{
            b1[j]=u1[j];
			b2[j]=u1[j]-c1[j];
			b3[j]=u1[j]+c1[j];
			}
  	        static double b1z[N+1],b1f[N+1],b2z[N+1],b2f[N+1],b3z[N+1],b3f[N+1];//定义分裂变量
		    for(j=0;j<N+1;j++)//为分裂变量赋值
			{
            b1z[j]=(b1[j]+sqrt(b1[j]*b1[j]+e*e))*1.0/2;
            b1f[j]=(b1[j]-sqrt(b1[j]*b1[j]+e*e))*1.0/2;
			b2z[j]=(b2[j]+sqrt(b2[j]*b2[j]+e*e))*1.0/2;
			b2f[j]=(b2[j]-sqrt(b2[j]*b2[j]+e*e))*1.0/2;
			b3z[j]=(b3[j]+sqrt(b3[j]*b3[j]+e*e))*1.0/2;
			b3f[j]=(b3[j]-sqrt(b3[j]*b3[j]+e*e))*1.0/2;
			}
		    static double fz1b[N+1],fz2b[N+1],fz3b[N+1];
		    static double ff1b[N+1],ff2b[N+1],ff3b[N+1];
		    static double wzb[N+1],wfb[N+1];
		    for(j=0;j<N+1;j++)//为每一个fz及ff赋值
			{
            fz1b[j]=d1[j]*1.0/2/r*(2*(r-1)*b1z[j]+b2z[j]+b3z[j]);
            ff1b[j]=d1[j]*1.0/2/r*(2*(r-1)*b1f[j]+b2f[j]+b3f[j]);
            fz2b[j]=d1[j]*1.0/2/r*(2*(r-1)*b1z[j]*u1[j]+b2z[j]*(u1[j]-c1[j])+b3z[j]*(u1[j]+c1[j]));
            ff2b[j]=d1[j]*1.0/2/r*(2*(r-1)*b1f[j]*u1[j]+b2f[j]*(u1[j]-c1[j])+b3f[j]*(u1[j]+c1[j]));
            wzb[j]=(3-r)*(b2z[j]+b3z[j])*c1[j]*c1[j]*1.0/2/(r-1);
            wfb[j]=(3-r)*(b2f[j]+b3f[j])*c1[j]*c1[j]*1.0/2/(r-1);
            fz3b[j]=d1[j]*1.0/2/r*((r-1)*b1z[j]*u1[j]*u1[j]+b2z[j]*(u1[j]-c1[j])*(u1[j]-c1[j])*1.0/2+b3z[j]*(u1[j]+c1[j])*(u1[j]+c1[j])*1.0/2+wzb[j]);
            ff3b[j]=d1[j]*1.0/2/r*((r-1)*b1f[j]*u1[j]*u1[j]+b2f[j]*(u1[j]-c1[j])*(u1[j]-c1[j])*1.0/2+b3f[j]*(u1[j]+c1[j])*(u1[j]+c1[j])*1.0/2+wfb[j]);
			}
			static double Lz1b[N+1],Lf1b[N+1],Lz2b[N+1],Lf2b[N+1],Lz3b[N+1],Lf3b[N+1];
			static double L1b[N+1],L2b[N+1],L3b[N+1];
			static double u2a[N+1],u2b[N+1],u2c[N+1];//这些分别对应R-K方法中的中间量u2
        fzd=FZ(fz1b);//调用WENO差分函数
		for(k=0;k<N+1;k++)
		{
			Lz1b[k]=0-fzd[k];          
		}
		ffd=FF(ff1b);
		for(k=0;k<N+1;k++)
		{
			Lf1b[k]=0-ffd[k];          
		}
		for(k=0;k<N+1;k++)
		{
			L1b[k]=Lz1b[k]+Lf1b[k];
		}
		fzd=FZ(fz2b);
		for(k=0;k<N+1;k++)
		{
			Lz2b[k]=0-fzd[k];          
		}
		ffd=FF(ff2b);
		for(k=0;k<N+1;k++)
		{
			Lf2b[k]=0-ffd[k];          
		}
		for(k=0;k<N+1;k++)
		{
			L2b[k]=Lz2b[k]+Lf2b[k];
		}
		fzd=FZ(fz3b);
		for(k=0;k<N+1;k++)
		{
			Lz3b[k]=0-fzd[k];          
		}
		ffd=FF(ff3b);
		for(k=0;k<N+1;k++)
		{
			Lf3b[k]=0-ffd[k];          
		}
		for(k=0;k<N+1;k++)
		{
			L3b[k]=Lz3b[k]+Lf3b[k];
		}
		for(k=0;k<N+1;k++)//第2个循环结构，为了求出数组u(2)[N+1]
		{
			u2a[k]=3.0/4*U1[k]+1.0/4*(u1a[k]+dt*L1b[k]);
			u2b[k]=3.0/4*U2[k]+1.0/4*(u1b[k]+dt*L2b[k]);
			u2c[k]=3.0/4*U3[k]+1.0/4*(u1c[k]+dt*L3b[k]);
		}//3阶R-K方法第二步结束
            static double u2[N+1],d2[N+1],p2[N+1];//分别对应R-K方法第一步产生的新的速度、密度、压力
		    for(j=0;j<N+1;j++)//根据u1a[k],u1b[k],u1c[k]产生新的u1[k],d1[k],p1[k];
			{
				d2[j]=u2a[j];
				u2[j]=u2b[j]*1.0/u2a[j];
				p2[j]=(u2c[j]-d2[j]*u2[j]*u2[j]*1.0/2)*(r-1);
			}
			static double c2[N+1];
			for(j=0;j<N+1;j++)//重新计算声速
			{
                c2[j]=sqrt(r*p2[j]*1.0/d2[j]);
			}
		    static double h1[N+1],h2[N+1],h3[N+1];//定义3个特征值
	        for(j=0;j<N+1;j++)//先通过一个循环结构产生N+1个特征值
			{
            h1[j]=u2[j];
			h2[j]=u2[j]-c2[j];
			h3[j]=u2[j]+c2[j];
			}
  	        static double h1z[N+1],h1f[N+1],h2z[N+1],h2f[N+1],h3z[N+1],h3f[N+1];//定义分裂变量
		    for(j=0;j<N+1;j++)//为分裂变量赋值
			{
            h1z[j]=(h1[j]+sqrt(h1[j]*h1[j]+e*e))*1.0/2;
            h1f[j]=(h1[j]-sqrt(h1[j]*h1[j]+e*e))*1.0/2;
			h2z[j]=(h2[j]+sqrt(h2[j]*h2[j]+e*e))*1.0/2;
			h2f[j]=(h2[j]-sqrt(h2[j]*h2[j]+e*e))*1.0/2;
			h3z[j]=(h3[j]+sqrt(h3[j]*h3[j]+e*e))*1.0/2;
			h3f[j]=(h3[j]-sqrt(h3[j]*h3[j]+e*e))*1.0/2;
			}
		    static double fz1h[N+1],fz2h[N+1],fz3h[N+1];
		    static double ff1h[N+1],ff2h[N+1],ff3h[N+1];
		    static double wzh[N+1],wfh[N+1];
		    for(j=0;j<N+1;j++)//为每一个fz及ff赋值
			{
            fz1h[j]=d2[j]*1.0/2/r*(2*(r-1)*h1z[j]+h2z[j]+h3z[j]);
            ff1h[j]=d2[j]*1.0/2/r*(2*(r-1)*h1f[j]+h2f[j]+h3f[j]);
            fz2h[j]=d2[j]*1.0/2/r*(2*(r-1)*h1z[j]*u2[j]+h2z[j]*(u2[j]-c2[j])+h3z[j]*(u2[j]+c2[j]));
            ff2h[j]=d2[j]*1.0/2/r*(2*(r-1)*h1f[j]*u2[j]+h2f[j]*(u2[j]-c2[j])+h3f[j]*(u2[j]+c2[j]));
            wzh[j]=(3-r)*(h2z[j]+h3z[j])*c2[j]*c2[j]*1.0/2/(r-1);
            wfh[j]=(3-r)*(h2f[j]+h3f[j])*c2[j]*c2[j]*1.0/2/(r-1);
            fz3h[j]=d2[j]*1.0/2/r*((r-1)*h1z[j]*u2[j]*u2[j]+h2z[j]*(u2[j]-c2[j])*(u2[j]-c2[j])*1.0/2+h3z[j]*(u2[j]+c2[j])*(u2[j]+c2[j])*1.0/2+wzh[j]);
            ff3h[j]=d2[j]*1.0/2/r*((r-1)*h1f[j]*u2[j]*u2[j]+h2f[j]*(u2[j]-c2[j])*(u2[j]-c2[j])*1.0/2+h3f[j]*(u2[j]+c2[j])*(u2[j]+c2[j])*1.0/2+wfh[j]);
			}
			static double Lz1h[N+1],Lf1h[N+1],Lz2h[N+1],Lf2h[N+1],Lz3h[N+1],Lf3h[N+1];
			static double L1h[N+1],L2h[N+1],L3h[N+1];
			//double u2a[N+1],u2b[N+1],u2c[N+1];//这些分别对应R-K方法中的中间量u2
        fzd=FZ(fz1h);//调用WENO差分函数
		for(k=0;k<N+1;k++)
		{
			Lz1h[k]=0-fzd[k];          
		}
		ffd=FF(ff1h);
		for(k=0;k<N+1;k++)
		{
			Lf1h[k]=0-ffd[k];          
		}
		for(k=0;k<N+1;k++)
		{
			L1h[k]=Lz1h[k]+Lf1h[k];
		}
		fzd=FZ(fz2h);
		for(k=0;k<N+1;k++)
		{
			Lz2h[k]=0-fzd[k];          
		}
		ffd=FF(ff2h);
		for(k=0;k<N+1;k++)
		{
			Lf2h[k]=0-ffd[k];          
		}
		for(k=0;k<N+1;k++)
		{
			L2h[k]=Lz2h[k]+Lf2h[k];
		}
		fzd=FZ(fz3h);
		for(k=0;k<N+1;k++)
		{
			Lz3h[k]=0-fzd[k];          
		}
		ffd=FF(ff3h);
		for(k=0;k<N+1;k++)
		{
			Lf3h[k]=0-ffd[k];          
		}
		for(k=0;k<N+1;k++)
		{
			L3h[k]=Lz3h[k]+Lf3h[k];
		}
		for(k=0;k<N+1;k++)//第2个循环结构，为了求出数组u2[N+1]
		{
			U1[k]=3.0/4*U1[k]+1.0/4*(u2a[k]+dt*L1h[k]);
			U2[k]=3.0/4*U2[k]+1.0/4*(u2b[k]+dt*L2h[k]);
			U3[k]=3.0/4*U3[k]+1.0/4*(u2c[k]+dt*L3h[k]);
		}//3阶R-K方法结束
		for(j=0;j<N+1;j++)
		{
			d[j]=U1[j];
			u[j]=U2[j]*1.0/U1[j];
			p[j]=(U3[j]-d[j]*u[j]*u[j]*1.0/2)*(r-1);
		}
	}//时间循环结束	
	    FILE *fp;//定义文件指针
		fp=fopen("c:WENOP.txt","w");
		int i=0;
		for(i=0;i<N+1;i++)//采用循环结构写文件
		{
	    fprintf(fp,"%20.10e%20.10e%20.10e%20.10e\n",i*1.0/N,fabs(u[i]),d[i],p[i]);
		}
		fclose(fp);
		cout<<"计算结果已保存在您设定的文件中"<<endl;
}
