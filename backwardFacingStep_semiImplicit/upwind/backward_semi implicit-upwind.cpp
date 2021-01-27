#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;
double dx,dy,dt,ust[2000][2000],un[2000][2000],un1[2000][2000],vst[2000][2000],vn[2000][2000],vn1[2000][2000],p[2000][2000];
double Vp,ru,rv,r2,Cu,Cv,F,uij,vij,bij,e,Re,rms,t,T,E;
int Nx,Ny;
double C1,C2,C3,C4;
double Fe,Fw,Fs,Fn;
void IC()
{

    for(int i=0;i<Nx;i++)
    {
        for(int j=0;j<Ny;j++)
        {
            un[i][j]=0.5;
            vn[i][j]=0.5;
            p[i][j]=0.1;
            ust[i][j]=0.5;
            vst[i][j]=0.5;
        }
    }

}

void velocity_convective1(int i,int j)
{
                F=0;
                Fe=((un[i+1][j]+un[i][j])/2)*dy;
                if(Fe>=0)
                {
                    C1=ust[i][j]*Fe;
                    F=Fe;
                }
                else
                {
                    C1=ust[i+1][j]*Fe;
                }
                 Fw=((un[i-1][j]+un[i][j])/2)*(-1*dy);
                if(Fw>=0)
                {
                    C2=ust[i][j]*Fw;
                    F=F+Fw;
                }
                else
                {
                    C2=ust[i-1][j]*Fw;
                }
                Fn=((vn[i][j+1]+vn[i][j])/2)*dx;
                if(Fn>=0)
                {
                    C3=ust[i][j]*Fn;
                    F=F+Fn;
                }
                else
                {
                    C3=ust[i][j+1]*Fn;
                    }
                    Fs=((vn[i][j-1]+vn[i][j])/2)*(-1*dx);
                if(Fs>=0)
                {
                    C4=ust[i][j]*Fs;
                    F=F+Fs;
                }
                else
                {
                    C4=ust[i][j-1]*Fs;
                }
                Cu=C1+C2+C3+C4;

}

void velocity_convective2(int i,int j)
{
                F=0;
                Fe=((un[i+1][j]+un[i][j])/2)*dy;
                if(Fe>=0)
                {
                    C1=vst[i][j]*Fe;
                    F=Fe;
                }
                else
                {
                    C1=vst[i+1][j]*Fe;
                }
                 Fw=((un[i-1][j]+un[i][j])/2)*(-1*dy);
                if(Fw>=0)
                {
                    C2=vst[i][j]*Fw;
                    F=F+Fw;
                }
                else
                {
                    C2=vst[i-1][j]*Fw;
                }
                Fn=((vn[i][j+1]+vn[i][j])/2)*dx;
                if(Fn>=0)
                {
                    C3=vst[i][j]*Fn;
                    F=F+Fn;
                }
                else
                {
                    C3=vst[i][j+1]*Fn;
                    }
                    Fs=((vn[i][j-1]+vn[i][j])/2)*(-1*dx);
                if(Fs>=0)
                {
                    C4=vst[i][j]*Fs;
                    F=F+Fs;
                }
                else
                {
                    C4=vst[i][j-1]*Fs;
                }
                Cv=C1+C2+C3+C4;

}


void Interior_velocity()
{
    for(int i=1;i<(Nx-1);i++)
        {
            for(int j=1;j<(Ny-1);j++)
            {
                velocity_convective1(i,j);

               ru=-(Vp/dt)*ust[i][j]+(Vp/dt)*un[i][j]+(Vp/Re)*((ust[i+1][j]-2*ust[i][j]+ust[i-1][j])/(dx*dx)+(ust[i][j+1]-2*ust[i][j]+ust[i][j-1])/(dy*dy))-Cu;
               uij=Vp/dt+(2/Re)*(1/(dx*dx)+1/(dy*dy))*Vp+F;
               ust[i][j]=ust[i][j]+ru/uij;
                velocity_convective2(i,j);

               rv=-(Vp/dt)*vst[i][j]+(Vp/dt)*vn[i][j]+(Vp/Re)*((vst[i+1][j]-2*vst[i][j]+vst[i-1][j])/(dx*dx)+(vst[i][j+1]-2*vst[i][j]+vst[i][j-1])/(dy*dy))-Cv;
               vij=Vp/dt+(2/Re)*(1/(dx*dx)+1/(dy*dy))*Vp+F;
               vst[i][j]=vst[i][j]+rv/vij;
               e=e+ru*ru+rv*rv;

            }
        }

}

void Boundary_velocity()
{
    //inlet
    double n=0.5;

    for(int j=(Ny/2);j<Ny;j++)
    {
        ust[0][j]=2*(36*n-24*n*n-12)-ust[1][j];
        n=n+dy;
        vst[0][j]=-vst[1][j];
    }

    //left boundary
    for (int j=0;j<(Ny/2);j++)
    {
        ust[0][j]=-ust[1][j];
        vst[0][j]=-vst[1][j];
    }
    //top boundary
    for(int i=1;i<(Nx-1);i++)
    {
        ust[i][Ny-1]=-ust[i][Ny-2];
        vst[i][Ny-1]=-vst[i][Ny-2];

    }
    //bottom boundary
    for(int i=1;i<(Nx-1);i++)
    {
        ust[i][0]=-ust[i][1];
        vst[i][0]=-vst[i][1];
    }


    //outlet boundary
    for(int j=0;j<Ny;j++)
    {
        ust[Nx-1][j]=ust[Nx-2][j];
        vst[Nx-1][j]=vst[Nx-2][j];
    }

}

void Interior_pressure()
{
    for(int i=1;i<(Nx-1);i++)
    {
        for(int j=1;j<(Ny-1);j++)
        {
            r2=(ust[i+1][j]-ust[i-1][j])/(2*dx)+(vst[i][j+1]-vst[i][j-1])/(2*dy)-dt*((p[i-1][j]-2*p[i][j]+p[i+1][j])/(dx*dx)+(p[i][j+1]-2*p[i][j]+p[i][j-1])/(dy*dy));
            bij=-2*dt*(1.0/(dx*dx)+1.0/(dy*dy));
            p[i][j]=p[i][j]+2*r2/bij;
            rms=rms+(r2)*(r2);

        }
    }
}

void  Boundary_pressure()
{
    //inlet

    for(int j=0;j<Ny;j++)
    {
        p[0][j]=p[1][j];
    }
//top
    for(int i=1;i<(Nx-1);i++)
    {
        p[i][Ny-1]=p[i][Ny-2];
    }
    //bottom
    for(int i=1;i<(Nx-1);i++)
    {
        p[i][0]=p[i][1];
    }
    //outlet
    for(int j=0;j<Ny;j++)
    {
        p[Nx-1][j]=-p[Nx-2][j];
    }

}
void update_velocity_BC()
{
        double n=0.5;
        for(int j=(Ny/2);j<Ny;j++)
    {
        un1[0][j]=2*(36*n-24*n*n-12)-un1[1][j];
        n=n+dy;
        vn1[0][j]=-vn1[1][j];
    }

    //left boundary
    for (int j=0;j<(Ny/2);j++)
    {
        un1[0][j]=-un1[1][j];
        vn1[0][j]=-vn1[1][j];
    }
    //top boundary
    for(int i=1;i<(Nx-1);i++)
    {
        un1[i][Ny-1]=-un1[i][Ny-2];
        vn1[i][Ny-1]=-vn1[i][Ny-2];

    }
    //bottom boundary
    for(int i=1;i<(Nx-1);i++)
    {
        un1[i][0]=-un1[i][1];
        vn1[i][0]=-vn1[i][1];
    }


    //outlet boundary
    for(int j=0;j<Ny;j++)
    {
        un1[Nx-1][j]=un1[Nx-2][j];
        vn1[Nx-1][j]=vn1[Nx-2][j];
    }

}

void update_velocity()
{
    for(int i=1;i<Nx-1;i++)
    {
        for(int j=1;j<Ny-1;j++)
        {
           un1[i][j]=ust[i][j]-dt*(p[i+1][j]-p[i-1][j])/(2*dx);
           vn1[i][j]=vst[i][j]-dt*(p[i][j+1]-p[i][j-1])/(2*dy);

        }
    }
    update_velocity_BC();

}


void transfer()
{
    T=0;
    for(int i=0;i<Nx;i++)
    {
        for(int j=0;j<Ny;j++)
            {

                T=T+((un1[i][j]-un[i][j])*(un1[i][j]-un[i][j])+(vn1[i][j]-vn[i][j])*(vn1[i][j]-vn[i][j]));
                un[i][j]=un1[i][j];
                vn[i][j]=vn1[i][j];

            }
    }


}


void display()
{
    ofstream out("backward semi explicit.txt");
    streambuf *coutbuf = std::cout.rdbuf(); //save old buf
    cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
    for(int j=0;j<Ny;j++)
    {
        for(int i=0;i<Nx;i++)
            {
               cout<<un1[i][j]<<'\t';

            }
            cout<<endl;
    }

    std::cout.rdbuf(coutbuf);


}

void TimeStep()
{
    while(T>0.001)
{
    e=1;
    while(e>0.000001)
    {   e=0;
        Interior_velocity();
        Boundary_velocity();
        e=sqrt(e/(2*Nx*Ny));
        //cout<<e<<endl;
}



    e=1;
    while(e>0.000001)
    {
        rms=0;
        Interior_pressure();
        Boundary_pressure();
        e=sqrt(rms/((Nx-2)*(Ny-2)));



    }
    update_velocity();
    transfer();

    T=sqrt(T/(2*Nx*Ny));
    T=T/dt;
    cout<<T<<endl;

}

    display();
}
void recirculation()
{
        double rel=0;

        for(int i=1;i<Nx;i++)
        {

            if(un1[i][1]<0)
            rel=rel+dx;
            else
            rel=rel+0;
        }
        cout<<rel;

}


int main()
{
    Nx=100.0;
    Ny=50.0;
    dx=10.0/(Nx-1);
    dy=1.0/(Ny-1);
    dt=0.01;
    Re=100;
    Vp=dx*dy;


    IC();
    T=1;


    TimeStep();

    recirculation();
    return 0;

}
