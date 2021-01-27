#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;
double dx,dy,dt,w1[2000][2000],w[2000][2000],u[2000][2000],v[2000][2000],si[2000][2000],si1[2000][2000];
double Vp,r1,r2,C,F,aij,bij,e,Re,rms,t,T,E;
int Nx,Ny;
double C1,C2,C3,C4;
double Fe,Fw,Fs,Fn;
void IC()
{

    for(int i=0;i<Nx;i++)
    {
        for(int j=0;j<Ny;j++)
        {
            w1[i][j]=0.5;
            w[i][j]=0.5;
        }
    }
    for(int i=0;i<Nx;i++)
    {
        for(int j=0;j<Ny;j++)
        {

            si[i][j]=0.5;
            si1[i][j]=0.5;
        }
    }





    for(int j=(Ny/2);j<Ny;j++)
    {
        u[0][j]=1.0;
        v[0][j]=0.0;

    }
    for(int i=1;i<Nx;i++)
    {
        for(int j=(Ny/2);j<Ny;j++)
        {
            u[i][j]=0.5;
            v[i][j]=0.0;

        }
    }
    for(int i=0;i<Nx;i++)
    {
        for(int j=0;j<(Ny/2);j++)
        {
            u[i][j]=0.5;
            v[i][j]=0.0;

        }
    }





}
void omega_convective(int i,int j)
{
                F=0;
                Fe=((u[i+1][j]+u[i][j])/2)*dy;
                if(Fe>=0)
                {
                    C1=w1[i][j]*Fe;
                    F=Fe;
                }
                else
                {
                    C1=w1[i+1][j]*Fe;
                    }
                 Fw=((u[i-1][j]+u[i][j])/2)*(-1*dy);
                if(Fw>=0)
                {
                    C2=w1[i][j]*Fw;
                    F=F+Fw;
                }
                else
                {
                    C2=w1[i-1][j]*Fw;
                }
                Fn=((v[i][j+1]+v[i][j])/2)*dx;
                if(Fn>=0)
                {
                    C3=w1[i][j]*Fn;
                    F=F+Fn;
                }
                else
                {
                    C3=w1[i][j+1]*Fn;
                    }
                    Fs=((v[i][j-1]+v[i][j])/2)*(-1*dx);
                if(Fs>=0)
                {
                    C4=w1[i][j]*Fs;
                    F=F+Fs;
                }
                else
                {
                    C4=w1[i][j-1]*Fs;
                }
                    C=C1+C2+C3+C4;
}


void Interior_omega()
{
    for(int i=1;i<(Nx-1);i++)
        {
            for(int j=1;j<(Ny-1);j++)
            {
                omega_convective(i,j);

                r1=Vp*(w[i][j])/dt-Vp*(w1[i][j]/dt-(1/Re)*((w1[i+1][j]-2*w1[i][j]+w1[i-1][j])/(dx*dx)+(w1[i][j+1]-2*w1[i][j]+w1[i][j-1])/(dy*dy)))-C;
                aij=Vp/dt + Vp/Re*(2/(dx*dx)+2/(dy*dy))+F;
                w1[i][j]=w1[i][j]+r1/aij;
                rms=rms+(r1)*(r1);

            }
        }
}

void Boundary_omega()
{
    //inlet
    double n=0.5;
    double c=0;
    for(int j=(Ny/2);j<Ny;j++)
    {
        c=-8*n*n*n+18*n*n-12*n+2.5;
        w1[0][j]=2*((8*(c-si[1][j]))/(dx*dx)+48*n-36)-w1[1][j];
        n=n+dy;
    }

    //left boundary
    for (int j=0;j<(Ny/2);j++)
    {
        w1[0][j]=-2*((8*si[1][j])/(dx*dx))-w1[1][j];
    }
    //top boundary
    for(int i=1;i<(Nx-1);i++)
    {
        w1[i][Ny-1]=2*(8*(0.5-si[i][Ny-2]))/(dy*dy)-w1[i][Ny-2];

    }
    //bottom boundary
    for(int i=1;i<(Nx-1);i++)
    {
        w1[i][0]=-2*((8*si[i][1])/(dy*dy))-w1[i][1];
    }


    //outlet boundary
    for(int j=0;j<Ny;j++)
    {
        w1[Nx-1][j]=w1[Nx-2][j];
    }

}

void Interior_psi()
{
    for(int i=1;i<(Nx-1);i++)
    {
        for(int j=1;j<(Ny-1);j++)
        {
            r2=-1*(w1[i][j]+(si1[i+1][j]-2*si1[i][j]+si1[i-1][j])/(dx*dx)+(si1[i][j+1]-2*si1[i][j]+si1[i][j-1])/(dy*dy));
            bij=-2*(1/(dx*dx)+1/(dy*dy));
            si1[i][j]=si1[i][j]+r2/bij;
            rms=rms+(r2)*(r2);


        }
    }
}

void  Boundary_psi()
{
    //inlet
    double n=0.5;
    double c=0;
    for(int j=(Ny/2);j<Ny;j++)
    {
        c=-8*n*n*n+18*n*n-12*n+2.5;
        si1[0][j]=2*c-si1[1][j];
        n=n+dy;
    }
    //left
    for(int j=0;j<(Ny/2);j++)
    {
        si1[0][j]=-si1[1][j];
    }
    //top
    for(int i=1;i<(Nx-1);i++)
    {
        si1[i][Ny-1]=1-si1[i][Ny-2];
    }
    //bottom
    for(int i=1;i<(Nx-1);i++)
    {
        si1[i][0]=-si1[i][1];
    }
    //outlet
    for(int j=0;j<Ny;j++)
    {
        si1[Nx-1][j]=si1[Nx-2][j];
    }

}

void velocity()
{
    for(int i=1;i<Nx-1;i++)
    {
        for(int j=1;j<Ny-1;j++)
            {
                u[i][j]=(si1[i][j+1]-si1[i][j-1])/(2*dy);
            }
    }
    for(int j=0;j<Ny/2;j++)
    {
            u[0][j]=-u[1][j];

    }
    for(int i=0;i<Nx;i++)
    {
        u[i][0]=-u[i][1];
        u[i][Ny-1]=-u[i][Ny-2];
    }
    for(int j=0;j<Ny;j++)
    {
        u[Nx-1][j]=u[Nx-2][j];
    }

   double n=0.5;
    for(int j=Ny/2;j<Ny;j++)
    {
        u[0][j]=2*(36*n-24*n*n-12)-u[1][j];
        n=n+dy;

    }


    for(int j=1;j<Ny-1;j++)
    {
        for(int i=1;i<(Nx-1);i++)
            {
                v[i][j]=-(si1[i+1][j]-si1[i-1][j])/(2*dx);
            }
    }
    for(int i=0;i<Nx;i++)
    {
            v[i][0]=-v[i][1];
            v[i][Ny-1]=-v[i][Ny-2];
    }
    for(int j=0;j<Ny;j++)
    {
        v[0][j]=-v[1][j];
        v[Nx-1][j]=v[Nx-2][j];
    }




}
void transfer()
{
    for(int i=0;i<Nx;i++)
    {
        for(int j=0;j<Ny;j++)
            {

                si[i][j]=si1[i][j];
                w[i][j]=w1[i][j];

            }
    }


}


void display()
{
    ofstream out("acfd_ass_1.txt");
    streambuf *coutbuf = std::cout.rdbuf(); //save old buf
    cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
    for(int j=0;j<Ny;j++)
    {
        for(int i=0;i<Nx;i++)
            {
               cout<<u[i][j]<<'\t';

            }
            cout<<endl;
    }
    std::cout.rdbuf(coutbuf);


}

void TimeStep()
{
    while(T<=10)
    {
        e=1;

    while(e>0.00001)
    {
        rms=0;
        Interior_omega();
        Boundary_omega();

        e=sqrt(rms/((Nx-2)*(Ny-2)));
       cout<<e;
      cout<<endl;

    }
    e=1;
    while(e>0.00001)
    {
        rms=0;
        Interior_psi();
        Boundary_psi();
        e=sqrt(rms/((Nx-2)*(Ny-2)));


    }
    velocity();
    transfer();
     T=T+dt;
}

    display();
}
void recirculation()
{
        double rel=0;

        for(int i=1;i<Nx;i++)
        {

            if(u[i][1]<0)
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
    Re=800;
    Vp=dx*dy;

    IC();
    T=dt;


    TimeStep();

    recirculation();
    return 0;

}
