//No sé muy bien por qué; pero el caso es que conseguí que ejecutara archivos de códigos en .cpp 
//Las reglas de oro para que funcione son 2: 
//          1) Que la carpeta donde tenga el código esté alojada en una carpeta dentro de la carpeta "usaurios" y no dentro de otra como "documentos" porque mac bloquea el acceso incluso habiendo reseteado todos los permisos de VsCode y habiendo especificado que puede usarlos
//          2) Cada vez que ejecute NO se guarda automáticamente el archivo así que manualmente debo guardar siempre
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

using namespace std;
//Defino constantes que usaré (En SI)
#define G 1.0
#define PI 3.14159265359
#define Mt 1.0//Masa de la Tierra
#define Ml 1.0//Masa de la Luna
#define Msol 1.0//Masa del Sol
#define Mmarte 1.0//Masa de marte

//Todas estas que incluyo aquí son las que debo usar según el problema que trate si es con tierra-luna o sistema solar
#define dTL 1.0//Distancia Tierra-Luna en afelio(La necesitaré para redimensionar y que los números no sean deamasiado diferentes y dé problemas numéricos)
#define dTS 1.0//Distancia Tierra-Sol en afelio (Para redimensionar)


//Es mejor definir las funciones de cada K para cada variable antes y luego sustituir los valores que necesites con Runge kutta (TODO LO QUE USO ESTÁ YA TILDADO)
double Calc_rpunto (vector<double> Pr, double h, int i);
double Calc_phipunto (vector<double> Pphi, vector<double> r, double h, int i);
double Calc_Prpunto (vector<double> r, vector<double> phi, vector<double> Pphi, vector<double> m, double h, int i);
double Calc_Pphipunto (vector<double> r, vector<double> phi, vector<double> m, double h, int i);

//Defino una funcion para sumar columnas de matriz a un vector
vector<double> Suma (vector<double> V, vector<vector<double> >M, int k);

//Quiero ver si se conserva el hamiltoniano
double H, H0, T, V, DeltaH;
double denV1 ,denV2, denV3, denT1, denT2, denT3;



int main()
{
    //Abro ficheros para almacenar datgos 
    ofstream M0, M1, M2, Hamiltoniano;

    M0.open("Masa0.txt");
    M1.open("Masa1.txt");
    M2.open("Masa2.txt");
    Hamiltoniano.open("Hamiltoniano.txt");

    if (M0.is_open()==true&&M1.is_open()==true&&M2.is_open()==true&&Hamiltoniano.is_open()==true)
    {
        //Comienzo definiendo todos los grados de libertad del problema (OJO TODAS LAS VARIABLES SON LAS TILDADAS AUNQUE NO LO INDIQUE)
        vector<double> r(3), phi(3), Pr(3), Pphi(3), velocidad_0x(3), velocidad_0y(3), x(3), y(3), x0(3), y0(3);//Necesitaré 12 ecuaciones diferenciales de primer orden
        vector<double> m(3);

        //Ahora declaro las variables y parámetros para hacer RK-4 (Cada K son matrices de forma que las filas son la el k 1, 2, 3 o 4 y las columnas la variable r, Pr o phi y Pphi)
        vector<vector<double> > Kr(4, vector<double>(3)), Kphi(4, vector<double>(3)), KPr(4, vector<double>(3)), KPphi(4, vector<double>(3));
        vector<double> Kr_arg (3), KPr_arg(3), Kphi_arg(3), KPphi_arg(3);
        //double h=0.00001;//Es el paso de tiempos en minutos PARA FIGURAS QUE NECESITEN MUCHA PRECISION COMO DRAGONFLY O BUTTERFLY
        double h=0.01; //Eso es para casos más estables como la figure 8. Para Broucke 0.001
        int i, j, k;//Contador para calcular nuevos valores de posiciones y tal

        
        //CONDICIONES INICIALES TIERRA-LUNA: se empezó a probar por tanto eran las distancia reales Tierra-Luna y la velocidad de la luna respecto a la horizontal en el afelio de su órbita
        //VELOCIDADES
        //velocidad_0y[0]=0;
        //velocidad_0y[1]=1000;
        //velocidad_0y[2]=-1000;
        //velocidad_0y[2]=30.8e3;//Velocidad luna si simulas junto a tierra
        
        //RADIOS
        //r[0]=0.001;
        //r[1]=dTS/dTS;
        //r[2]=dTS/dTS;
        //r[2]=(dTL+dTS)/dTS;


        //Para las trayectorias cerradas raras hay que liarla parda con las coordenadas porque con la restricción de que la velocidad inicial sea perpendicular al eje x no llegamos a nada
        //Declaro antes las coordenadas iniciales en coordenadas cartesianas para que sea fácil colocar cosas. ESTAS SON LAS  QUE DEBO CAMBIAR EN CADA ÓRBITA
        
        //1) PARA LA EIGHT-FIGURE CLASSIC ORBITA
        //COORDENADAS DE LAS POSICIONES. OJO QUE DESCENTRO PARA EVITAR SINGULARIDADES Y NO DIVIDIR POR 0
        x[0]=-5.0;
        y[0]=-5.0;

        x[1]=-5.0+0.97000436;
        y[1]=-5.0-0.24308753;

        x[2]=-5.0-0.97000436;
        y[2]=-5.0+0.24308753;

        //VELOCIDADES CARTESIANAS (ESTO ES LO QUE DEBO CAMBIAR EN CADA ORBITA)
        velocidad_0x[0]=-0.93240737;
        velocidad_0y[0]=-0.86473146;

        velocidad_0x[1]=-0.5*velocidad_0x[0];
        velocidad_0y[1]=-0.5*velocidad_0y[0];

        velocidad_0x[2]=velocidad_0x[1];
        velocidad_0y[2]=velocidad_0y[1];
        //////////////////////////////////

        //2) PARA LA FIGURE 8 V.1.B
        //COORDENADAS DE LAS POSICIONES
        //x[0] = 5+0.000001; // Un desfase infinitesimal para evitar r=0 exacto en el cuerpo central
        //y[0] = 5+0.0;

        //x[1] = 5-1.0;
        //y[1] = 5+0.0;

        //x[2] = 5+1.0;
        //y[2] = 5+0.0;

        // Valores numéricos puros de la tabla del paper
        //double v1 = 0.339393;
        //double v2 = 0.536191;

        // Cuerpo 0 (Central en el origen)
        //velocidad_0x[0] = -2.0 * v1;
        //velocidad_0y[0] = -2.0 * v2;

        // Cuerpo 1 (Compensamos el signo que le va a meter el cos(PI) y sin(PI) en la proyección)
        //velocidad_0x[1] = v1; 
        //velocidad_0y[1] = v2;

        // Cuerpo 2 (Se queda igual porque cos(0) = 1)
        //velocidad_0x[2] = v1;
        //velocidad_0y[2] = v2;
        //////////////////////////////////

        //3) PARA LA FIGURE 8 V.8.B
        //COORDENADAS DE LAS POSICIONES
        //x[0] = 5+0.000001; // Un desfase infinitesimal para evitar r=0 exacto en el cuerpo central
        //y[0] = 5+0.0;

        //x[1] = 5-1.0;
        //y[1] = 5+0.0;

        //x[2] = 5+1.0;
        //y[2] = 5+0.0;

        // Valores numéricos puros de la tabla del paper
        //double v1 = 0.151718;
        //double v2 = 0.5298;

        // Cuerpo 0 (Central en el origen)
        //velocidad_0x[0] = -2.0 * v1;
        //velocidad_0y[0] = -2.0 * v2;

        // Cuerpo 1 (Compensamos el signo que le va a meter el cos(PI) y sin(PI) en la proyección)
        //velocidad_0x[1] = v1; 
        //velocidad_0y[1] = v2;

        // Cuerpo 2 (Se queda igual porque cos(0) = 1)
        //velocidad_0x[2] = v1;
        //velocidad_0y[2] = v2;
        //////////////////////////////////

        //4) PARA LA DRAGONFLY 4.A
        //COORDENADAS DE LAS POSICIONES
        //x[0] = 5+0.000001; // Un desfase infinitesimal para evitar r=0 exacto en el cuerpo central
        //y[0] = 5+0.0;

        //x[1] = 5-1.0;
        //y[1] = 5+0.0;

        //x[2] = 5+1.0;
        //y[2] = 5+0.0;

        // Valores numéricos puros de la tabla del paper
        //double v1 = 0.080584;
        //double v2 = 0.588836;

        // Cuerpo 0 (Central en el origen)
        //velocidad_0x[0] = -2.0 * v1;
        //velocidad_0y[0] = -2.0 * v2;

        // Cuerpo 1 (Compensamos el signo que le va a meter el cos(PI) y sin(PI) en la proyección)
        //velocidad_0x[1] = v1; 
        //velocidad_0y[1] = v2;

        // Cuerpo 2 (Se queda igual porque cos(0) = 1)
        //velocidad_0x[2] = v1;
        //velocidad_0y[2] = v2;
        //////////////////////////////////

        //5) PARA LA BUTTERFLY IVB.3.A
        //COORDENADAS DE LAS POSICIONES
        //x[0] = 5+0.000001; // Un desfase infinitesimal para evitar r=0 exacto en el cuerpo central
        //y[0] = 5+0.0;

        //x[1] = 5-1.0;
        //y[1] = 5+0.0;

        //x[2] = 5+1.0;
        //y[2] = 5+0.0;

        // Valores numéricos puros de la tabla del paper
        //double v1 = 0.405916;
        //double v2 = 0.230163;

        // Cuerpo 0 (Central en el origen)
        //velocidad_0x[0] = -2.0 * v1;
        //velocidad_0y[0] = -2.0 * v2;

        // Cuerpo 1 (Compensamos el signo que le va a meter el cos(PI) y sin(PI) en la proyección)
        //velocidad_0x[1] = v1; 
        //velocidad_0y[1] = v2;

        // Cuerpo 2 (Se queda igual porque cos(0) = 1)
        //velocidad_0x[2] = v1;
        //velocidad_0y[2] = v2;
        //////////////////////////////////

         //6) PARA LA BROUCKE A2
        //COORDENADAS DE LAS POSICIONES
        //x[0]=0.336130095;
        //y[0]=0.0;

        //x[1]=0.7699893804;
        //y[1]=0.0;

        //x[2]=-1.1061194753;
        //y[2]=0.0;

        //VELOCIDADES CARTESIANAS (ESTO ES LO QUE DEBO CAMBIAR EN CADA ORBITA)
        //velocidad_0x[0]=0.0;
        //velocidad_0y[0]=1.532431537;

        //velocidad_0x[1]=0.0;
        //velocidad_0y[1]=-0.6287350978;

        //velocidad_0x[2]=0.0;
        //velocidad_0y[2]=-0.9036964391;
        //////////////////////////////////
       
        
        //Calculo los radios y ángulos iniciales en polares
        for (i=0; i<3; i++)
        {
            r[i]=sqrt(x[i]*x[i]+y[i]*y[i]);
            phi[i]=atan2(y[i],x[i]);
            cout<<"radios iniciales="<<r[i]<<" y ángulos iniciales="<<phi[i]*(180/PI)<<endl;
        }


        //Si desarrollas los momentos en r y phi con las derivadas de la velocidad generalizada del lagrangiano y sacas la expresión del Pr y Pphi, si haces después el reescalado y sustituyes r punto y phi punto usando la igualdad de V en cartesianas y en polares pero en funcion de los versores de cartesianas sale esto:
        for(i=0; i<3; i++)
        {
            Pr[i]=(velocidad_0y[i]*sin(phi[i])+velocidad_0x[i]*cos(phi[i]))/(dTL);
            Pphi[i]=(r[i]*(velocidad_0y[i]*cos(phi[i])-velocidad_0x[i]*sin(phi[i])))/(dTL*dTL);
        }
        
       
        m[0]=1.0;
        m[1]=1.0; 
        m[2]=1.0;



        //Almaceno el primer valor de la x y de y en los ficheros correspondientes
        M0<<x[0]<<'\t'<<y[0]<<endl;
        M1<<x[1]<<'\t'<<y[1]<<endl;
        M2<<x[2]<<'\t'<<y[2]<<endl;
        ///////////////////////////////////////////////////////////////////////////

        //El algoritmo de runge Kutta
        for (i=0; i<7000000; i++)//Bucle para el paso temporal
        {
            //Primero necesito que el K1 reciba de argumento las variables del paso anterior
            Kr_arg=r;
            Kphi_arg=phi;
            KPr_arg=Pr;
            KPphi_arg=Pphi;
            for (k=0; k<4; k++)//Bucle para recorrer los K desde 1 a 4
            {
                for (j=0; j<3; j++)//Bucle para recorrer los 3 cuerpos
                {
                    Kr[k][j]=Calc_rpunto(KPr_arg, h, j);
                    Kphi[k][j]=Calc_phipunto(KPphi_arg, Kr_arg, h, j);
                    KPr[k][j]=Calc_Prpunto(Kr_arg, Kphi_arg, KPphi_arg, m, h, j);
                    KPphi[k][j]=Calc_Pphipunto(Kr_arg, Kphi_arg, m, h, j);
                }
                //Ahora cambio los argumentos para que se les sume como toque el K que acabo de calcular
                Kr_arg=Suma(r, Kr, k);
                Kphi_arg=Suma(phi, Kphi, k);
                KPr_arg=Suma(Pr, KPr, k);
                KPphi_arg=Suma(Pphi, KPphi, k);
            }

            //Aquí debo calcular la nueva r, Pr, phi, Pphi. Hago un bucle que recorra cada cuerpo
            for(j=0; j<3; j++)
            {
                r[j]=r[j]+1.0/6.0*(Kr[0][j]+2*Kr[1][j]+2*Kr[2][j]+Kr[3][j]);
                Pr[j]=Pr[j]+1.0/6.0*(KPr[0][j]+2*KPr[1][j]+2*KPr[2][j]+KPr[3][j]);
                Pphi[j]=Pphi[j]+1.0/6.0*(KPphi[0][j]+2*KPphi[1][j]+2*KPphi[2][j]+KPphi[3][j]);
                phi[j]=phi[j]+1.0/6.0*(Kphi[0][j]+2*Kphi[1][j]+2*Kphi[2][j]+Kphi[3][j]);
            }
                
            //Ahora es sólo guardar los cálculos en variables x e y que es lo que realmente almaceno en ficheros
            for(j=0; j<3; j++)
            {
                x[j]=r[j]*cos(phi[j]);
                y[j]=r[j]*sin(phi[j]);
            }
            M0<<x[0]<<'\t'<<y[0]<<endl;
            M1<<x[1]<<'\t'<<y[1]<<endl;
            M2<<x[2]<<'\t'<<y[2]<<endl;

            //Lo que implemento ahora es para calcular el hamiltoniano y ver si es constante H=T+V

            //para calcular la energía cinética y que no diverja cuando pasen cerca del origen hago que se anulen los denominadores más pequeños de la cuenta
            denT1=r[0]*r[0];
            //if(denT1<1e-6)
            //{
                //denT1=1e6;
            //}

            denT2=r[1]*r[1];
            //if(denT2<1e-6)
            //{
                //denT2=1e6;
            //}

            denT3=r[2]*r[2];
            //if(denT3<1e-6)
            //{
                //denT3=1e6;
            //}

            T=((Pphi[0]*Pphi[0]/denT1+Pr[0]*Pr[0])*m[0]+(Pphi[1]*Pphi[1]/denT2+Pr[1]*Pr[1])*m[1]+(Pphi[2]*Pphi[2]/denT3+Pr[2]*Pr[2])*m[2])*(dTS*dTS/2.0);

            //para calcular la energía potencial en los denominadores me aseguro de que no salgan numeros demasiado chicos que disparen la energía
            denV1=sqrt(r[0]*r[0]-2*r[0]*r[1]*cos(phi[0]-phi[1])+r[1]*r[1]);
            //if(denV1<1e-6)
            //{
                //denV1=1e6;
            //}
            
            denV2=sqrt(r[0]*r[0]-2*r[0]*r[2]*cos(phi[0]-phi[2])+r[2]*r[2]);
            //if(denV2<1e-6)
            //{
                //denV2=1e6;
            //}
            
            denV3=sqrt(r[1]*r[1]-2*r[1]*r[2]*cos(phi[1]-phi[2])+r[2]*r[2]);
            //if(denV3<1e-6)
            //{
                //denV3=1e6;
            //}


            V=-G/(dTL)*((m[0]*m[1]/denV1)+(m[0]*m[2]/denV2)+(m[1]*m[2]/denV3));
            H=T+V;

            if(i==0)
            {
                H0=H;
            }

            DeltaH=H-H0;
            //if (abs(DeltaH/H0)>0.2)
            //{
                //DeltaH=0;
            //}
            Hamiltoniano<<H<<'\t'<<DeltaH/H0<<endl;
        }


    }

    else
    {
        cerr<<"El fichero no se abre (Gracias Bernier por enseñarnos la buena praxis qué haríamos sin ti)";
    }
    
    cout<<"Se ejecutó para Hamilton bien";
    
    return 0;
}

//Defino todas las funciones que uso después

double Calc_rpunto (vector<double> Pr, double h, int i)
{
    double solucion;
    
    solucion=Pr[i]*h;
    
    return solucion;
}

double Calc_phipunto (vector<double> Pphi, vector<double> r, double h, int i)
{
    
    if (r[i] < 1e-6) 
    {
        return 0.0;
    } 
    
    double solucion;
 
    solucion=Pphi[i]/(r[i]*r[i])*h;

    

    return solucion;
}

double Calc_Prpunto (vector<double> r, vector<double> phi, vector<double> Pphi, vector<double> m, double h, int i)//El contador i es para saber en qué cuerpo estoy y aplicar un cálculo u otro
{
    double Term[4], suma, den1, den2;

    int j = (i + 1) % 3; //Si i=0 entonces j=1. Si i=1 entonces j=2. Si i=2 entonces j=0
    int k = (i + 2) % 3; //Si i=0 entonces k=2. Si i=1 entonces k=0. Si i=2 entonces k=1

    

    den1=pow(r[i]*r[i] - 2*r[i]*r[j]*cos(phi[i]-phi[j]) + r[j]*r[j], 1.5);
    den2=pow(r[i]*r[i] - 2*r[i]*r[k]*cos(phi[i]-phi[k]) + r[k]*r[k], 1.5);

    if (den1<3.5e-7)
    {
        den1=1e8;
    }

    if(den2<3.5e-7)
    {
        den2=1e8;
    }

    Term[0] = -(G * r[i] * m[j]) / den1;
    Term[1] = -(G * r[i] * m[k]) / den2;
    Term[2] =  (G * r[j] * m[j] * cos(phi[i]-phi[j])) / den1;
    Term[3] =  (G * r[k] * m[k] * cos(phi[i]-phi[k])) / den2;

    //si r[i] es cero, el término Pphi^2 / r^3 no se añade
    double termino_centrifugo = 0.0;
    if (r[i] > 1e-6) 
    {
        termino_centrifugo = (Pphi[i] * Pphi[i]) / (r[i] * r[i] * r[i]);
    }

    suma = Term[0] + Term[1] + Term[2] + Term[3] + termino_centrifugo;

    return suma * h;
}


double Calc_Pphipunto (vector<double> r, vector<double> phi, vector<double> m, double h, int i)
{
    double Term[2], suma, den1, den2;

    int j = (i + 1) % 3;
    int k = (i + 2) % 3;

    den1=pow(r[i]*r[i] - 2*r[i]*r[j]*cos(phi[i]-phi[j]) + r[j]*r[j], 1.5);
    den2=pow(r[i]*r[i] - 2*r[i]*r[k]*cos(phi[i]-phi[k]) + r[k]*r[k], 1.5);

    if (den1<3.5e-7)
    {
        den1=1e8;
    }
    
    if (den2<3.5e-7)
    {
        den2=1e8;
    }

    //El torque ejercido por los otros dos cuerpos sobre la masa i
    Term[0] = (r[j] * m[j] * sin(phi[i] - phi[j])) / den1;
    Term[1] = (r[k] * m[k] * sin(phi[i] - phi[k])) / den2;

    suma = Term[0] + Term[1];
    
    return -G * r[i] * suma * h;
}

vector<double> Suma (vector<double>V, vector<vector<double> >M, int k)//LA funcion va a sumar al array de 3 dimensiones de cada variable la mitad o el K correspondiente del paso anterior. Como K es una matriz, necesito una funcion que sume a un vector las columnas de una matriz. A ESTA FUNCION SE IGUALA LOS Ki_arg DEL MAIN
{
    int l;
    vector<double> SUMA(3);
    
    if (k==0 || k==1)//En este caso sumamos la mitad del K que toque
    {
        for (l=0; l<3; l++)
        {
            SUMA[l]=V[l]+1.0/2.0*M[k][l];
        }
    }
    else if (k==2)
    {
        for (l=0; l<3; l++)
        {   
            SUMA[l]=V[l]+M[k][l];
        }
    }
    return SUMA;
}