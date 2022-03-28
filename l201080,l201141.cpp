//RUNTIME TERROR
// A Project by: Nabeeha Mudassir and Nisa Nadeem
#include<iostream>
#include<fstream>
#include<string>
#include<cstdlib>
#include<cmath>

using namespace std;

#define MaxRows 500
#define MaxCols 500

struct grayImage{
    grayImage(){
        Rows = Cols = 0;
        Loaded = 0;
        for(int r = 0; r< MaxRows; r++)
            for(int c = 0; c< MaxCols; c++)
                Image[r][c] = 0;
    }

    unsigned short setPixel(unsigned short value, int r, int c){
        if( r >= Rows || c >= Cols || r < 0 || c < 0)
            return -1;
        Image[r][c] = value;
        return value;
    }

    int getPixel(int r, int c){
        if( r >= Rows || c >= Cols || r < 0 || c < 0)
            return -1;
        return Image[r][c];
    }

    int setRows(int rows){
        if(rows < 1 || rows > MaxRows)
            return -1;
        Rows = rows;
        return Rows;
    }

    int getRows(){
        return Rows;
    }

    int setCols(int cols){
        if(cols < 1 || cols > MaxCols)
            return -1;
        Cols = cols;
        return Cols;
    }

    int getCols(){
        return Cols;
    }
    void fillcolor (int startrow, int endrow, int topcol, int bottomcol, int fillval ){
            for (int i=startrow;i<endrow;i++)
                {
                    for (int j=topcol;j<bottomcol;j++)
                        {
                            Image[i][j]=fillval;
                        }
                }
    }
    void combineSidebySide(grayImage &Two, int fillValue){
        if (Rows>Two.Rows)
            {       int startrow=Two.Rows;//where the second image ends
                    int endrow=Rows+Two.Rows;//boundary of fill area
                    int topcol=Cols;
                    int bottomcol=Cols+Two.Cols;//boundary of fill area

                    fillcolor(startrow,endrow,topcol,bottomcol,fillValue);
            }
        else if (Rows<Two.Rows)
        {
                int startrow=Rows;
                int endrow=Two.Rows;
                int topcol=0;
                int bottomcol=Cols;
                fillcolor(startrow,endrow,topcol,bottomcol,fillValue);
        }
        if (Rows>Two.Rows || Rows<Two.Rows || Rows==Two.Rows)
            {
                for (int i=0;i<Two.Rows;i++)
                    {
                        for (int j=0;j<Two.Cols;j++)
                        {
                           if (Cols+j<MaxCols)
                            {
                                Image[i][j+Cols]=Two.Image[i][j];
                            }
                        }
                    }
                Cols=Cols+Two.Cols;
                if (Cols>MaxCols)
                    {Cols=MaxCols;}
                if (Rows<Two.Rows)
                    {Rows=Two.Rows;}
            }
}
    void combineTopToBottom(grayImage &Two, int fillValue){
                if (Cols>Two.Cols)
                {
                    int startrow=Rows;
                    int endrow=Rows+Two.Rows;
                    int startcol=0;
                    int endcol=Cols+Two.Cols;
                    fillcolor(startrow,endrow,startcol,endcol,fillValue);
                }
                else if (Cols<Two.Cols)
                {
                    int startrow=0;
                    int endrow=Rows;
                    int startcol=Cols;
                    int endcol=Two.Cols;
                    fillcolor (startrow,endrow,startcol,endcol,fillValue);
                }
                if (Cols<Two.Cols || Cols>Two.Cols ||Cols==Two.Cols)
                        {
                                for (int i=0;i<Two.Rows;i++)
                                {
                                    for (int j=0;j<Two.Cols;j++)
                                    {
                                        if (Rows+i<MaxRows)
                                        {Image[i+Rows][j]=Two.Image[i][j];}
                                    }
                                }
                                Rows=Rows+Two.Rows;
                                if (Rows>MaxRows)
                                {
                                    Rows=MaxRows;
                                }
                                if (Cols<Two.Cols)
                                {
                                    Cols=Two.Cols;
                                }
                        }
}
    void fadein (grayImage &Two, int frames, int seconds, string BaseFileName){
            grayImage Fil;
            double step=1.0/(frames*seconds);
            int R=Rows;
            if (R<Two.Rows)
                {
                    R=Two.Rows;
                }
            int C=Cols;
            if (C<Two.Cols)
                {
                    C=Two.Cols;
                }

            Fil.Rows=R;
            Fil.Cols=C;
            Fil.Loaded=1;
            Fil.Maximum=Maximum;

            if (Maximum<Two.Maximum)
                {
                    Fil.Maximum=Two.Maximum;
                }
            int i=0;
            for (double avg=1;avg>=0;avg=avg-step)
            {
                    for (int row=0;row<R;row++)
                    {
                        for (int column=0;column<C;column++)
                        {
                            Fil.Image[row][column]=avg*Image[row][column] + (1-avg)*Two.Image[row][column];
                        }
                    }
                    char New[10];
                    itoa (i,New,10);//it converts integer to ascii when it is given: value to be converted, the new array to store,base that is ten for deci
                    Fil.Save(BaseFileName + New + ".pgm");
                    i=i+1;
                    if (avg-step < 0 && avg-step > -step)
                    {
                        avg=0;
                    }
            }
}
    void flip (int f){
        short unsigned int flipped[Rows][Cols];
        if (f==1)//1 is vertical
        {
            for (int r=0;r<Rows;r++)
            {
                for (int col=0;col<Cols;col++)
                {
                    flipped[r][col]=Image[r][Cols-1-col];
                }
            }

            for (int r=0;r<Rows;r++)
            {
                for (int col=0;col<Cols;col++)
                {
                    Image[r][col]=flipped[r][col];
                }
            }
        }
        if (f==0)//0 is horizontal
        {
            for (int r=0;r<Rows;r++)
            {
                for (int col=0;col<Cols;col++)
                {
                    flipped[Rows-1-r][col]=Image[r][col];
                }
            }

            for (int r=0;r<Rows;r++)
            {
                for (int col=0;col<Cols;col++)
                {
                    Image[r][col]=flipped[r][col];
                }
            }

        }

}
    void Rotate(double angle, int aboutx, int abouty, grayImage& RotatedImage){
        int inew,jnew;
        angle=(3.141592*angle)/180.0;//sine only works with radians
        for (int r=0;r<Rows;r++)
        {
            for (int c=0;c<Cols;c++)
            {
                inew=((r-aboutx)*cos(angle)-(c-abouty)*sin(angle) + aboutx);
                jnew=((r-aboutx)*sin(angle)+(c-abouty)*cos(angle)+ abouty);

                if (inew>0 && inew<MaxRows && jnew>0 && jnew<MaxCols)
                {
                    RotatedImage.Image[inew][jnew]=Image[r][c];
                }
            }
        }
        RotatedImage.Rows=MaxRows;
        RotatedImage.Cols=MaxCols;
        RotatedImage.Maximum=Maximum;
        RotatedImage.Loaded=1; //to satisfy the condition inside the save function
}
    void changeBrightness(int amount){
        if (amount<Maximum && amount>=0)
        {
            for (int i=0;i<Rows;i++)
            {
                for (int j=0;j<Cols;j++)
                {
                    if (Image[i][j]+amount<Maximum)
                    {Image[i][j]=Image[i][j]+amount;}

                    else {Image[i][j]=Maximum;}
                }
            }
        }
    }
    void Negative(){
        for (int i=0;i<Rows;i++)
        {
            for (int j=0;j<Cols;j++)
            {
                Image[i][j]=Maximum-Image[i][j];
            }
        }
    }
    void sortarray (int window[]){
        int temp, i , j;
        for(i = 0; i < 9; i++)
            {
                temp = window[i];
                for(j = i-1; j >= 0 && temp < window[j]; j--)
                {
                    window[j+1] = window[j];
                }
            window[j+1] = temp;
            }
}
    void medianFilter(grayImage &Result, int filterSize = 3){
       int n=filterSize*filterSize;
       int window[n]; //filter size can only be an odd number
       for (int i=1;i<Rows-1;i++)
        {
            for (int j=1;j<Cols-1;j++)
            {

                window[0]=Image[i-1][j-1];
                window[1]=Image[i][j-1];
                window[2]=Image[i+1][j-1];
                window[3]=Image[i-1][j];
                window[4]=Image[i][j];
                window[5]=Image[i+1][j];
                window[6]=Image[i-1][j+1];
                window[7]=Image[i][j+1];
                window[8]=Image[i+1][j+1];

                sortarray(window);
                Result.Image[i][j]=window[4];
                Result.Loaded=1;
            }
        }
}
    void meanFilter(grayImage& Result, double filterSize = 3){
        int window[9];
        for (int i=1;i<Rows-1;i++)
        {
            for (int j=1;j<Cols-1;j++)
            {
                window[0]=Image[i-1][j-1];
                window[1]=Image[i][j-1];
                window[2]=Image[i+1][j-1];
                window[3]=Image[i-1][j];
                window[4]=Image[i][j];
                window[5]=Image[i+1][j];
                window[6]=Image[i-1][j+1];
                window[7]=Image[i][j+1];
                window[8]=Image[i+1][j+1];

                int mean=(window[0]+window[1]+window[2]+window[3]+window[4]+window[5]+window[6]+window[7]+window[8])/9;
                Result.Image[i][j]=mean;
                Result.Loaded=1;
            }
        }
}
    void Filter(grayImage& Result,double Mask[3][3]){
        for (int r=1;r<Rows-1;r++)
        {
            for (int c=1;c<Cols-1;c++)
            {
                double sum=0.0;
                for(int k = -1; k <= 1;k++)
                    {
                        for(int j = -1; j <=1; j++)
                        {
                            sum = sum + Mask[j+1][k+1] * Image[r - j][c - k];
                        }

                    }
                    Result.Image[r][c]=sum;
            }
        }
    Result.Maximum=Maximum;
    Result.Cols=Cols;
    Result.Rows=Rows;
    Result.Loaded=1;

}
    void DerivativeImage(grayImage& R, grayImage& New) {
    double M[3][3]={  {-1,0,1}, {-1,0,1},  {-1,0,1} };//DOES NOT WORK WHY
    double N[3][3]={  {-1,-1,-1}, {0,0,0}, {1,1,1} };

    Filter(R,M);
    Filter(New,N);

    for (int i=1;i<Rows-1;i++)
    {
        for (int j=1;j<Cols-1;j++)
        {
            double a=(R.Image[i][j]*R.Image[i][j]) + (New.Image[i][j] * New.Image[i][j]);
            if (a<0)
                {
                    a=a*-1;
                }

            short unsigned  b=sqrt(a);

            if (b>=0 && b<=255)
                {
                    Image[i][j]= b;
                }
        }
    }
}
    void Crop (grayImage& Result, const int startrow,  const int endrow,  const int startcol, const int endcol, bool resizeornot= false ){
            if (startrow<0 || startrow>Rows || endrow<0 || endrow>Rows || startcol<0 || startcol>Cols || startcol<0 || endcol>Cols)
            {
                cout<<endl<<"you have entered something out of bound "<<endl;
                return;
            }

            int i=0;
            for ( int r=startrow ; r<=endrow ; r++ , i++)
            {
                int j=0;
                for (int c=startcol ; c<=endcol ; c++ , j++)
                {
                    Result.Image[i][j]=Image[r][c];
                }
            }
            Result.Cols=endcol-startcol;
            Result.Rows=endrow-startrow;
            Result.Maximum=Maximum;
            Result.Loaded=1;

}
    void Resize(grayImage& Result,double NewRows, double NewColumns){
        double skiprows;
        double skipcols;
        double Ratio=NewColumns/Cols;
       // cout<<"Ratio: "<<Ratio;
        if (Ratio<=1) //result image will be smaller
            {
                skiprows=1.0/Ratio;
                skipcols=1.0/Ratio;
                int R,C,oldc,oldr;
                double newr=0;
                for (double r=0  ; r<Rows ; newr++, r=r+skiprows)
                {
                    double newc=0;
                    for (double c=0 ; c<Cols ; c=c+skipcols , newc++)
                    {
                         R=newr;
                         C=newc;
                         oldc=c;
                         oldr=r;
                         Result.Image[R][C]=Image[oldr][oldc];
                    }
                }
            }
       else if (Ratio>1)
            {
               double sr=NewRows;
               sr/=Rows;
               double sc=NewColumns;
               sc/= Cols;
               int r=sr;
               int c=sc;

               for (int i=0;i<NewRows;i++)
               {
                   for (int j=0;j<NewColumns;j++)
                   {
                       int I,J;
                       I=i/sr;
                       J=j/sc;


                       Result.Image[i][j]=Image[I][J];
                   }
                   /*Result.Rows=NewRows;
                   Result.Cols=NewColumns;
                   Result.Maximum=Maximum;
                   Result.Loaded=1;*/
               }


            }
        Result.Loaded=1;
        Result.Cols=NewColumns;
        Result.Rows=NewRows;
        Result.Maximum=Maximum;
}
    void Resize(grayImage& Result, double Ratio){
        double newrows=Ratio*Rows;
        double newcols=Ratio*Cols;
        Resize(Result,newrows,newcols);
    }
    void Quantize (grayImage& Result, int Levels){
        if (Levels>0 && Levels<Maximum)
        {
            int number=Maximum/Levels;

            for (int r=0;r<Rows;r++)
            {
                for (int c=0;c<Cols;c++)
                {
                 Result.Image[r][c] = Image[r][c] / number * number +  number/2;
                }
            }
            Result.Loaded=1;
            Result.Maximum=Maximum;
            Result.Rows=Rows;
            Result.Cols=Cols;
        }
}
    void Transform(grayImage& Result , double Matrix[3][3]){

        for (int i=0;i<Rows;i++)
        {
            for (int j=0;j<Cols;j++)
            {
                double sumi,sumj,sumk;
                double newi,newj;
                int NI,NJ;
                sumi = i * Matrix[0][0] + j * Matrix[0][1] + Matrix[0][2];
                sumj = i * Matrix[1][0] + j * Matrix[1][1] + Matrix[1][2];
                sumk = i * Matrix[2][0] + j * Matrix[2][1] + Matrix[2][2];
                if (sumk!=0)
                { newi=sumi/sumk;
                  newj=sumj/sumk;
                  NI=newi;
                  NJ=newj;
                }
                if (newi<MaxRows && newj<MaxCols)
                Result.Image[NI][NJ]=Image[i][j];
            }
        }
        Result.Loaded=1;
        Result.Rows=Rows;
        Result.Cols=Cols;
        Result.Maximum=Maximum;
}

    int Save(string File_Name){

       if (Loaded!=1)
           {
               cout<<"The file was not loaded.";
               return 2;
           }
       else
        {ofstream Output(File_Name.c_str());
        Output<<"P2\n"<<"#RUNTIME TERROR: A Project by Nabeeha Mudassir and Nisa Nadeem "<<endl<<Cols<<" "<<Rows<<endl<<Maximum<<endl;
        for (int i=0;i<Rows;i++)
        {
            for (int j=0;j<Cols;j++)
            {
                Output<<Image[i][j]<<" ";
            }
            Output<<endl;
        }
        //cout<<endl<<"Saved Successfully."<<endl;
        Output.close();
        }

    }
    int load(string File_Name){
        if (File_Name.substr(File_Name.find_last_of(".")+1)!="pgm")
            {
               cout<<"file format not supported. Only pgm is allowed. ";
               return 1;
            }
        char MagicNum[10];
        char Comment[200];
        int cols,rows,MaxValue;

        ifstream File(File_Name.c_str());

        if (!File)
        {
            cout<<"Could not open file for reading! Please ensure this file actually exists. ";
            return 2;
        }

        else
        {File.getline(MagicNum,10);
        File.getline(Comment,100);
        File>>cols>>rows>>MaxValue;

        Maximum=MaxValue;
        Rows=rows;
        Cols=cols;

        setRows(rows);
        setCols(cols);

        for (int r=0;r<rows;r++)
        {
            for (int c=0;c<cols;c++)
            {
               File>>Image[r][c];
            }
        }
        /*cout<<MagicNum;
        cout<<endl<<Comment<<endl<<Cols<<" "<<Rows<<endl<<Maximum;*/
        Loaded=1;
        cout<<endl<<"Loaded Successfully."<<endl;
        File.close();
        }
}
private:
    unsigned short Image[MaxRows][MaxCols];
    int Rows, Cols, Maximum;
    int Loaded=0;
};

int main(){
    grayImage GM1,GM2,GM3;
    cout<<"RUNTIME TERROR: A Project by Nabeeha Mudassir and Nisa Nadeem"<<endl<<endl;
    cout<<endl<<"What would you like to do with your image? Choose an option from our menu. "<<endl<<endl;
    cout<<"Press 1 to LOAD an image."<<endl<< "Press 2 to SAVE this image."<<endl<<"Press 3 to FLIP it."<<endl<<"Press -3 to ROTATE it."<<endl;
    cout<<"Press 4 to change its BRIGHTNESS"<<endl<<"Press 5 to apply MEDIAN FILTER."<<endl<<"Press 6 to apply MEAN FILTER."<<endl;
    cout<<"Press 7 to apply simple FILTER"<<endl<<"Press 8 to find DERIVATIVE. "<<endl;
    cout<<"Press -8 to find NEGATIVE of an image."<<endl<<"Press 9 to CROP it."<<endl<<"Press 10 to RESIZE it."<<endl;
    cout<<"Press 11 to QUANTIZE it"<<endl<<"Press 12 to TRANSFORM it."<<endl;
    cout<<"Press 13 to apply FADE IN to 2 images."<<endl<<"Press 14 to combine these 2 images SIDE BY SIDE."<<endl;
    cout<<"Press 15 to combine 2 images TOP TO BOTTOM."<<endl;
    cout<<"If you want to exit the loop, then press 0."<<endl;

    int n;
    cout<<"Enter your choice: ";
    cin>>n;
    while (n!=0)
    {

        if ( n==1)
        {
            cout<<"Enter the name of the image, ensure .pgm is included. "<<endl;
            string image;
            cin>>image;
            GM1.load(image);
        }
        else if (n==2)
        {
            string img;
            cout<<"Name your saved image. Include .pgm extension."<<endl;
            cin>>img;
            GM1.Save(img);
        }
        else if (n == 3)
        {
            cout<<endl<<"Press 1 for vertical flip or 0 for horizontal flip."<<endl;
            int f;
            cin>>f;
            if ( f==1 || f==0 )
            {
                GM1.flip(f);
                cout<<endl<<"Flip complete. New image will be saved as flipped."<<endl;
                GM1.Save("flipped.pgm");
            }
            else cout<<endl<<"Invalid choice. "<<endl;
        }
        else if (n == -3)
        {
            double angle; int x; int y;
            cout<<endl<<"Enter the angle in degrees by which you would like to rotate."<<endl;
            cin>>angle;
            cout<<endl<<"Now Enter the x coordinate followed by y coordinate about which you would like to rotate."<<endl;
            cin>>x>>y;
            cout<<endl<<"Rotated file will be saved in another file called Rotated."<<endl;
            GM1.Rotate(angle,x,y,GM2);
            GM2.Save("RotatedFile.pgm");
        }
        else if (n == 4)
        {
            cout<<endl<<"Enter the amount by which you would like to change the brightness."<<endl;
            int amount;
            cin>>amount;
            GM1.changeBrightness(amount);
            cout<<endl<<"New file will be saved as Brighter"<<endl;
            GM1.Save("brighterr.pgm");
        }
        else if (n == 5)
        {
            GM1.medianFilter(GM1,3);
            cout<<endl<<"New image is saved as MedianF"<<endl;
            GM1.Save("MedianF.pgm");
        }
        else if (n == 6)
        {
            GM1.meanFilter(GM1,3);
            cout<<endl<<"New image is saved as MeanF"<<endl;
            GM1.Save("MeanF.pgm");
        }
        else if (n == 7)
        {
            cout<<endl<<"You will need to enter a 3 by 3 mask. Enter 9 values. "<<endl;
            double Mask [3][3];
            for (int i=0;i<3;i++)
            {
                for (int j=0;j<3;j++)
                {
                    cin>>Mask[i][j];
                }
            }
            GM1.Filter(GM2,Mask);
            cout<<endl<<"New File is saved as Filtered."<<endl;
            GM2.Save("filtered.pgm");
        }
        else if (n == 8)
        {
            GM1.DerivativeImage(GM2,GM3);
            cout<<endl<<"New File will be saved as Derived"<<endl;
            GM1.Save("Derived.pgm");
        }
        else if (n == -8)
        {
            GM1.Negative();
            cout<<endl<<"New image will be saved as Negativer"<<endl;
            GM1.Save("Negativer.pgm");
        }
        else if (n == 9)
        {
            int r1,r2,c1,c2,y;
            cout<<endl<<"Enter the dimensions you want in your cropped image."<<endl;
            cout<<"Enter the start row, followed by ending row, then enter starting column, followed by ending column."<<endl;
            cin>>r1>>r2>>c1>>c2;
            cout<<"If you would like to Resize the Cropped Image as well, then press 1, else press 0."<<endl;
            cin>>y;
            if (y == 1)
            {
                GM1.Crop(GM2,r1,r2,c1,c2,true);
                int r, c;
                cout<<endl<<"Enter the new rows, followed by the new columns in the resized image."<<endl;
                cin>>r>>c;
                GM2.Resize(GM2,r,c);
                cout<<endl<<"New image will be saved as KroppedandResized"<<endl;
                GM2.Save("KroppedandResized.pgm");
            }
            else if (y == 0)
            {
                GM1.Crop(GM2,r1,r2,c1,c2,false);
                cout<<endl<<"Cropped image will be saved as Kropped."<<endl;
                GM2.Save("Kropped.pgm");
            }

        }
        else if (n == 10)
        {
                int r;
                cout<<endl<<"Press 1 to resize by ratio or press 2 to Resize by NewRows and NewCols"<<endl;
                cin>>r;
                if (r==1)
                {   double ratioo;
                    cout<<endl<<"Enter Ratio"<<endl;
                    cin>>ratioo;
                    GM1.Resize(GM2,ratioo);
                    cout<<endl<<"New Image will be saved as ResizedbyRatio"<<endl;
                    GM2.Save("ResizedbyRatio.pgm");
                }
                else if (r==2)
                 {
                    int r, c;
                    cout<<endl<<"Enter the new rows, followed by the new columns in the resized image."<<endl;
                    cin>>r>>c;
                    GM1.Resize(GM2,r,c);
                    cout<<endl<<"New image will be saved as Resized."<<endl;
                    GM2.Save("Resized.pgm");
                 }
                 else cout<<endl<<"Only enter 1 or 2."<<endl;
        }
        else if ( n== 11)
        {
            int l;
            cout<<endl<<"Enter the levels you want to quantize in"<<endl;
            cin>>l;
            GM1.Quantize(GM2,l);
            cout<<"New image will be saved as Quantized"<<endl;
            GM1.Save("Quantized.pgm");
        }
        else if (n == 12)
        {
            double Matrix[3][3];
            cout<<endl<<"Enter the 3 by 3 transformation matrix"<<endl;
            for (int i=0;i<3;i++)
            {
                for (int j=0;j<3;j++)
                    cin>>Matrix[i][j];
            }
            GM1.Transform(GM2,Matrix);
            cout<<endl<<"New image will be saved as Transformed."<<endl;
            GM1.Save("Transformed.pgm");
        }
        else if (n == 13)
        {
            double frame, seconds;
            string f,s;
            cout<<endl<<"Enter the first image name"<<endl;
            cin>>f;
            GM1.load(f);
            cout<<endl<<"Enter the second image name"<<endl;
            cin>>s;
            GM2.load(s);
            cout<<endl<<"Enter frames followed by seconds."<<endl;
            cin>>frame>>seconds;
            cout<<endl<<"Create a folder called FADED in your directory beforehand. So images can be saved in it."<<endl;
            GM1.fadein(GM2,frame,seconds,"Faded\\Images\\");
        }
        else if (n == 14)
        {
            string first;
            string second;
            cout<<endl<<"Enter the first image name. "<<endl;
            cin>>first;
            GM1.load(first);
            cout<<endl<<"Enter the second image name. "<<endl;
            cin>>second;
            GM2.load(second);
            cout<<endl<<"If images are not the same size, remaining space will be filled by a color. Choose a color. Enter a value from 0 to 255."<<endl;
            int fillc;
            cin>>fillc;
            GM1.combineSidebySide(GM2,fillc);
            cout<<endl<<"New image will be saved as SidebySide"<<endl;
            GM1.Save("SidebySide.pgm");
        }
        else if (n == 15)
        {
            string f,s;
            cout<<endl<<"Enter the first image name"<<endl;
            cin>>f;
            GM1.load(f);
            cout<<endl<<"Enter the second image name"<<endl;
            cin>>s;
            GM2.load(s);
            cout<<endl<<"If images are not the same size, remaining space will be filled by a color. Choose a color. Enter a value from 0 to 255."<<endl;
            int fillc;
            cin>>fillc;
            GM1.combineTopToBottom(GM2,fillc);
            cout<<endl<<"New image will be saved as ToptoBottom."<<endl;
            GM1.Save("ToptoBottom.pgm");
        }
        else
            {cout<<endl<<"Only enter from the choices specified."<<endl;}

    cout<<endl<<"Your choice: ";
    cin>>n;
    }
    if (n == 0)
        cout<<endl<<"You have terminated the loop. Program will end now.";

    return 0;

}
