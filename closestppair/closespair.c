#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <dos.h>
#include <math.h>

#include <float.h>

#define M_PI 3.14159265358979323846
#define DIM 3

typedef struct  {

	double coordinates[DIM];
	double coordinates_projected[DIM];
	double dist_projection_to_origine;

} Point;

typedef struct  {
	Point p;
	Point q;
	double min_distance;

} Result;

int counting=0;
Point *p_array;

Point *A[DIM];

Point *arr1;
Point *arr2;
char *paths[200];


char return_path_[200];

void array_copy(Point x[], Point y[],int size_, int begin);

void fill_point_array() ;
void fill_projected_point_array() ;
void fill_AA();

void print_array(Point x[]);
void print_array_with_size(Point x[],int size__);



void merge_sort(Point arr[],int low,int high,int sort_type_flag,int coordinate_no);
void merge(Point arr[],int l,int m,int h,int sort_type_flag,int coordinate_no);

Result smartForce(Point *A_[],int flag);
Result bruteForceClosestPair(Point p_array_[],int size_,int type);
Result recursiveClosestPair(Point A_[],int _size);

void random_points();

void recursivin_smartforcu_gectigi_nokta(int limit,int dim,int scope,int flag);
void Histogram(int flag_);
void Histogram_projected(int flag_);

void print_result(Result r);


FILE *myfile;
double *AA;

double calculate_distance(Point p,Point q);

double * convert_to_xyz(double *coordinates);
void convert_to_uvw(int i);

double basis[DIM][DIM];
double *inv[DIM]={0};
double origine[DIM];
double detrm( double a[DIM][DIM], double k );
double** cofact( double num[ DIM ][ DIM ], double f );
double** trans( double num[ DIM ][ DIM ], double fac[ DIM ][ DIM ], double r );
void generate_basis();
double* projection_vector(double *u, double *v);
double *sum_of_vectors(double *u, double *v);
double *substruction_of_vectors(double *u, double *v);
int flag__=0;
int flag2__=0;

static int size=100;
int  recursivin_smartforcu_gectigi_nokta_=2;
int change_size(int i);
 int a=0;
int main(void)
{ //int N;
	double t1,t2; double total_time=0; int count=0; double total1=0;double total2=0;double total3=0; int i;
	int size_user=0;  double determinant; int choice=0; int change_size_value=1;

	Result r1;
	Result r2;
	Result r3;
	printf("\nDimension: %d\n",DIM);
	while(change_size_value==1 )
	{

	    //change_size(change_size_value);

		p_array= (Point *) malloc(sizeof(Point)* size);
		arr1= (Point *) malloc(sizeof(Point)* size);
		arr2= (Point *) malloc(sizeof(Point)* size);
		AA = (double *) malloc(sizeof(Point) * size*DIM) ;


		choice=2;
		printf("-----------------------------------------------\n");
		switch(choice)
		{
		case 1:
			strcpy ( return_path_, "ornek.txt" );
			fill_point_array() ;
			break;
		case 2:
			strcpy ( return_path_, "ornek.txt" );
			random_points();
			break;
		case 3:

            // strcpy ( return_path_, "C:\\Users\\PBell\\Documents\\Visual Studio 2008\\Projects\\closestPair_k_dim\\closestPair_k_dim\\noktalar_db\\short.txt" );
			 strcpy ( return_path_, "C:\\Users\\PBell\\Documents\\Visual Studio 2008\\Projects\\closestPair_k_dim\\closestPair_k_dim\\long.txt");

//C:\Users\PBell\Documents\Visual Studio 2008\Projects\closestPair_k_dim\closestPair_k_dim\noktalar_db
			fill_AA();
			break;


		}



		/*
		for(i=0;i<DIM;i++)
		{
			if(i%2==0)
		{
		    recursivin_smartforcu_gectigi_nokta(recursivin_smartforcu_gectigi_nokta_,i,size/2,1);
		}
		else
		{
		    recursivin_smartforcu_gectigi_nokta(recursivin_smartforcu_gectigi_nokta_,i,size/2,3);
		}

		}
		*/



		for(count=0; count<1; count++)
		{

			t2=0; t1=0;

			//-----------------------------------------------

			//if(change_size_value<3)
			{
			t1=(double) clock();
			//printf ("It took me %d clicks (%f seconds).\n",t1,(double)t1/CLOCKS_PER_SEC);
			r1=bruteForceClosestPair(p_array,size,0);
			t2 =(double) clock();
			//printf ("It took me %d clicks (%f seconds).\n",t2,(double)t2/CLOCKS_PER_SEC);
			total1=total1+ t2-t1;
			}
			//else
			{
			//if(flag__==0) printf("Bruteforce runtime cok uzun\n\n");
			flag__=1;
			}




			//printf ("bruteForceClosestPair: %.f clicks (%f seconds).\n\n",t2-t1,((double)t2-(double)t1)/CLOCKS_PER_SEC);

			t2=0; t1=0;

			//-----------------------------------------------
			t1=(double) clock();

//for(i=0;i<DIM;i++)
//origine[i]=0;


//generate_basis();

            /*
			for(i=0;i<DIM;i++)
			inv[i]= (double *) malloc(sizeof(double)* DIM);

			 determinant = detrm( basis, DIM );
            // printf( "THE DETERMINANT IS=%f", determinant );

             if ( determinant == 0 )
             {//  printf( "\nMATRIX IS NOT INVERSIBLE\n" );
			 }
              else
             cofact( basis, DIM );
			  */

//fill_projected_point_array() ;


			//printf("p_array: \n\n");
			//print_array(p_array);





			for(i=0;i<DIM;i++)
			{
				A[i]= (Point *) malloc(sizeof(Point)* size);
				memcpy ( A[i], p_array,sizeof(Point)* size );

			}

			//for(i=0;i<DIM;i++)
			{
			//Histogram( i);
			//Histogram_projected(i);
			}
			flag2__=1;

			for(i=0;i<DIM;i++)
			{
				merge_sort(A[i],0,size-1, 1,i);
				//printf("A[%d]: \n\n",i);
				//print_array(A[i]);
			}

			//print_array(p_array);



			r2=smartForce(A,1);
			t2 =(double) clock();
			total2=total2+ t2-t1;

			for(i=0;i<DIM;i++)
			{
				free(A[i]);

			}

			//printf ("It took me %d clicks (%f seconds).\n",t2,(double)t2/CLOCKS_PER_SEC);

			//printf ("smartForce: %f clicks (%f seconds).\n\n",t2-t1,((double)t2-(double)t1)/CLOCKS_PER_SEC);


			t2=0; t1=0;

			//-----------------------------------------------


	t1= (double)clock();
			for(i=0;i<DIM;i++)
			{
				A[i]= (Point *) malloc(sizeof(Point)* size);
				memcpy ( A[i], p_array,sizeof(Point)* size );

				merge_sort(A[i],0,size-1, 1,i);
				//printf("A[%d]: \n\n",i);
				//print_array(A[i]);
			}



			r3=recursiveClosestPair(A[0] , size);
			t2 =(double) clock();
			total3=total3+ t2-t1;

			for(i=0;i<DIM;i++)
			{
				free(A[i]);

			}


			//printf ("recursiveClosestPair: %f clicks (%f seconds).\n",t2-t1,((double)t2-(double)t1)/CLOCKS_PER_SEC);

			//printf("\n\n\n------------------------------------------\n\n\n");
			t2=0; t1=0;


			//print_array(p_array	);



		}
		    flag2__=0;

			//if(flag__==0)
			{
			printf("bruteForceClosestPair'da en yakin iki double nokta: ");
			print_result(r1);
			printf ("bruteForceClosestPair: %f clicks (%f seconds).\n\n",total1/count,total1/(CLOCKS_PER_SEC*count));

			}


			printf("quickcp'de en yakin iki double nokta: ");
			print_result(r2);
			printf ("quickcp: %f clicks (%f seconds).\n\n",total2/count,total2/(CLOCKS_PER_SEC*count));

			printf("recursiveClosestPair'da en yakin iki double nokta: ");
			print_result(r3);
			printf ("recursiveClosestPair: %f clicks (%f seconds).\n\n",total3/count,total3/(CLOCKS_PER_SEC*count));




		free(p_array);
		free(arr1);
		free(arr2);

		total1=0;
		total2=0;
		total3=0;


		change_size_value++;
		}

	scanf("%d",&a);

	return 0;
}

int change_size(int i)
{
	switch (i)
	{
	    case 1: size=2000;
		break;

		case 2: size=3000;
		break;

		case 3: size=5100;
		break;

	    case 4: size=15000;
		break;

		case 5: size=85000;
		break;

		case 6: size=150000;
		break;


	}

}

void fill_AA()
{

	int i=0;
	int index_=0; int m=0; int j=0; int k=0;double temp=0;

    char line[50];

	//printf("%s\n",return_path_);

	if((myfile = fopen(return_path_, "r"))==NULL)
	{
	printf("Cannot open fileeeee.\n");

	}


	while( i<size*DIM)
	{
	fgets(line, 50, myfile) ;
	// printf("%s\n",line);
	 temp =(double) atof(line);
	if(temp!=0 && temp <1000 && temp>-1000)
	{
		AA[i]= (double) atof(line);
	//printf("AA[%d]: %.15f\n",i,AA[i]);
	i++;

	}


	}

	printf("Number of points: %d\n\n",size);
	fclose(myfile);   //close the file prior to exiting the routine

	i=0;


	for(j=0;j<size*DIM;j++)
	{

		for(m=0;m<DIM;m++)
		{
			p_array[k].coordinates[m]=AA[j+m];
			 //printf("%.15f ",p_array[k].coordinates[m]);
		}


		// printf("\n");


		k++;
		j=j+DIM-1;
	}


//	printf("number of distinct doubles: %d\n",j);
//	print_array(p_array);


	free(AA);


}


double *sum_of_vectors(double *u, double *v)
{
  double *projected_vector={0}; int i=0;
  projected_vector=(double *) malloc(sizeof(double)*DIM);

  for(i=0;i<DIM;i++)
	  projected_vector[i]=u[i]+v[i];

  return projected_vector;
}


double *substruction_of_vectors(double *u, double *v)
{
  double *projected_vector={0}; int i=0;
  projected_vector=(double *) malloc(sizeof(double)*DIM);

  for(i=0;i<DIM;i++)
	  projected_vector[i]=u[i]-v[i];

  return projected_vector;
}




double* projection_vector(double *u, double *v)
{
	double *projected_vector={0};int i=0;
	double inner_product=0; double norm_u=0;
    projected_vector=(double *) malloc(sizeof(double)*DIM);

	for(i=0;i<DIM;i++)
	{
	  inner_product=inner_product+u[i] *v[i];
	  norm_u=norm_u+ u[i]*u[i];
	}
	for(i=0;i<DIM;i++)
	{
	  projected_vector[i]= u[i] * (inner_product/norm_u);
	}

	return projected_vector;
}
void generate_basis()
{
	int i=0;int j=0;int k=0;int l=0; double norm=0; double u1,u2; int temp_int;
	double standard_basis[DIM][DIM]; double temp_basis[DIM][DIM];double temp_vector[DIM];double coefficients[DIM];

     double *temp_;

	 srand((unsigned int)time(NULL));
	 temp_=(double *) malloc(sizeof(double)*DIM);
	for(i=0;i<DIM;i++)
	{
	  coefficients[i]=(double)rand()/(M_PI *100);
	  temp_vector[i]=0;
	  temp_[i]=0;
	}

	for(i=0;i<DIM;i++)
	{
	   for(j=0;j<DIM;j++)
	   {
	      if(i==j)
		  {
		    standard_basis[i][j]=1;
			temp_basis[i][j]=0;
		  }
		  else
		  {
		    standard_basis[i][j]=0;
			temp_basis[i][j]=0;
		  }


	   }


	}

	if(DIM==2)
	{
	u1= (double)rand()/(M_PI *100);
	temp_int=u1;
	u1=u1-temp_int;
	basis[0][0]=u1;
	u2=1-u1*u1;
	u2=sqrtf(u2);
	basis[0][1]=u2;


	basis[1][0]=-u2;

	basis[1][1]=u1;


	}
	else{
	for(i=0;i<DIM;i++)
	{
		for(k=0;k<DIM;k++)
		{

			for(j=0;j<DIM;j++)
			{
				if(i==j)
				{
					continue;
				}
				else
				{
					temp_vector[k]= temp_vector[k]+standard_basis[j][k]*coefficients[j];

				}


			}
		}



		for(l=0;l<DIM;l++)
		{
			temp_basis[i][l]=temp_vector[l];
			temp_vector[l]=0;
		}


	}

	for(k=0;k<DIM;k++)
	basis[0][k]=temp_basis[0][k];

	for(k=1;k<DIM;k++)
	{
	  for(j=0;j<k;j++)
	  {
	    temp_=sum_of_vectors(temp_, projection_vector(basis[j],temp_basis[k]));
	  }


	  for(i=0;i<DIM;i++)
	  {
	   basis[k][i]=temp_basis[k][i]-temp_[i] ;

	  }

	  for(l=0;l<DIM;l++)
		  temp_[l]=0;

	}


	for(i=0;i<DIM;i++)
	{
		coefficients[i]=0;
	  for(j=0;j<DIM;j++)
	  {
	     coefficients[i]=coefficients[i]+ basis[i][j] * basis[i][j];

	  }
	    coefficients[i] =sqrtf(coefficients[i]);
	}

	for(i=0;i<DIM;i++)
	{
	  for(j=0;j<DIM;j++)
	  {
	    basis[i][j] = basis[i][j]/coefficients[i];
	  }
	}

	}


	  printf("\nORTHONORMAL BASIS: \n");

	for(i=0;i<DIM;i++)
	{
	  for(j=0;j<DIM;j++)
	  {
	    printf("%.4f ",basis[i][j]);

	  }
	  printf("\n");
	}
	  printf("\n");

}

double * convert_to_xyz(double *coordinates)
{
	double *converted_coordinates={0};
	int i=0; int j=0;
    converted_coordinates=(double *) malloc(sizeof(double)*DIM);
	for(i=0;i<DIM;i++)
	{

		converted_coordinates[i]=0;
	}

	//for(i=0;i<DIM;i++)
	//	printf("\n%.4f \n",coordinates[i]);
	for(i=0;i<DIM;i++)
	{
	   for(j=0;j<DIM;j++)
	   {
	    converted_coordinates[i]=converted_coordinates[i]+basis[i][j]*(coordinates[j]-origine[j]);

	   }

	}

return converted_coordinates;
}
void convert_to_uvw(int i)
{

int index=0; int j=0;

	for(index=0;index<DIM;index++)
	{
		p_array[i].coordinates_projected[index]=0;
	}

	for(index=0;index<DIM;index++)
	{
	   for(j=0;j<DIM;j++)
	   {
	    p_array[i].coordinates_projected[index]=p_array[i].coordinates_projected[index]+basis[index][j]*(p_array[i].coordinates[j]-origine[j]);

	   }

	}



}


double detrm( double a[DIM][DIM], double k )
{
    double s = 1, det = 0, b[DIM][DIM];
    int i, j, m, n, c;

    if ( k == 1 )
        {
            return ( a[ 0 ][ 0 ] );
        }
    else
        {
            det = 0;

            for ( c = 0;c < k;c++ )
                {
                    m = 0;
                    n = 0;

                    for ( i = 0;i < k;i++ )
                        {
                            for ( j = 0;j < k;j++ )
                                {
                                    b[ i ][ j ] = 0;

                                    if ( i != 0 && j != c )
                                        {
                                            b[ m ][ n ] = a[ i ][ j ];

                                            if ( n < k-2 )
											 n++;
                                            else
                                                {
                                                    n = 0;
                                                    m++;
                                                }
                                        }
                                }
                        }

                    det = det + s * ( a[ 0 ][ c ] * detrm( b, k-1 ) );
                    s = -1 * s;
                }
        }

    return ( det );
}

double** cofact( double num[ DIM ][ DIM ], double f )
{
    double b[ DIM ][ DIM ], fac[ DIM ][ DIM ]; int summ;
    int p, q, m, n, i, j;

    for ( q = 0;q < f;q++ )
        {
            for ( p = 0;p < f;p++ )
                {
                    m = 0;
                    n = 0;

                    for ( i = 0;i < f;i++ )
                        {
                            for ( j = 0;j < f;j++ )
                                {
                                    b[ i ][ j ] = 0;

                                    if ( i != q && j != p )
                                        {
                                            b[ m ][ n ] = num[ i ][ j ];

                                            if ( n < f-2 )
                                                n++;
                                            else
                                                {
                                                    n = 0;
                                                    m++;
                                                }
                                        }
                                }
                        }
                    summ=q+p;
                    fac[ q ][ p ] = pow((long double)(-1), summ ) * detrm( b, f-1 );
                }
        }

   return trans( num, fac, f );
}

double** trans( double num[ DIM ][ DIM ], double fac[ DIM ][ DIM ], double r )

{
    int i, j;
    double b[ DIM ][ DIM ], d;



    for ( i = 0;i < r;i++ )
        {
            for ( j = 0;j < r;j++ )
                {
                    b[ i ][ j ] = fac[ j ][ i ];
                }
        }

    d = detrm( num, r );
   // inv[ i ][ j ] = 0;

    for ( i = 0;i < r;i++ )
        {
            for ( j = 0;j < r;j++ )
                {
                    inv[ i ][ j ] = b[ i ][ j ] / d;
                }
        }

	 printf( "\nTHE MATRIX:\n" );

    for ( i = 0;i < r;i++ )
        {
            for ( j = 0;j < r;j++ )
                {
                    printf( "\t%f", basis[ i ][ j ] );
                }

            printf( "\n" );
        }


    printf( "\nTHE INVERSE OF THE MATRIX:\n" );

    for ( i = 0;i < r;i++ )
        {
            for ( j = 0;j < r;j++ )
                {
                    printf( "\t%f", inv[ i ][ j ] );
                }

            printf( "\n" );
        }


	return inv;
}


void print_result(Result r)
{
	int i=0;

	printf(" ( ");
	for(i=0;i<DIM;i++)
		printf(" %.15f ",  r.p.coordinates[i]);

	printf(") and (");


	for(i=0;i<DIM;i++)
		printf(" %.15f ",  r.q.coordinates[i]);
	printf(" )  ");

	printf("\naradaki uzaklik: %.15f\n", r.min_distance);




}



void recursivin_smartforcu_gectigi_nokta(int limit,int dim,int scope, int flag)
{

	int i=0; int j=1;
	if(flag==1)
	{
	for(i=0;i<scope;i++)
	{
		     p_array[i].coordinates[dim]=(double) j;

		    j++;


	if(j==limit)
	j=1;
	}


	}

	if(flag==3)
	{

	for(i=size-1;i>=scope;i--)
	{
		     p_array[i].coordinates[dim]=(double) j;


		    j++;

	if(j==limit)
	j=1;
	}


	}

}


void array_copy(Point *x, Point *y,int size_, int begin)
{

	int i=0; int n= size_-begin; int sum=0;
	for(i=0;i<n  ;i++)
	{
		sum=i+begin;
		x[i]=y[sum];
		//printf("i: %d sum:%d \n",i,sum);
	}
}
void print_array_with_size(Point x[],int size__)
{



	int i=0; int j=0;

	while(j<size__)
	{
		printf("( ");
		for(i=0;i<DIM;i++)
			printf(" %.4f ",  x[j].coordinates[i]);
		printf(") \n");
		/*
		printf("projected:( ");
		for(i=0;i<DIM;i++)
			printf(" %.4f ",  x[j].coordinates_projected[i]);
		printf(") \n");

		printf("( ");
		for(i=0;i<DIM;i++)
			printf("projected_vertical: %.4f ",  x[j].coordinates_projected_vertical[i]);
		printf(") \n");
		*/
		j++;
	}


	printf("\n");

}

void print_array(Point x[])
{



	int i=0; int j=0;

	while(j<size)
	{
		printf("( ");
		for(i=0;i<DIM;i++)
			printf(" %.8f ",  x[j].coordinates[i]);
		printf(") ");

		printf("projected: ( ");
		for(i=0;i<DIM;i++)
			printf(" %.8f ",  x[j].coordinates_projected[i]);
		printf(") \n");


		j++;
	}


	printf("\n");
}




void fill_point_array()
{

	int i=0;
	int index_=0; int m=0;
	char * double_string;
	char line[150];

	//printf("%s\n",return_path_);
	if((myfile = fopen(return_path_, "r"))==NULL)
	{
	printf("Cannot open fileeeee.\n");

	}


	while( i<size)
	{
	fgets(line, 150, myfile) ;
	//// printf("%s\n",line);


	double_string=strtok (line," ");
	p_array[i].coordinates[m] = atof(double_string);
	//printf("i: %d = %.6f ",i, p_array[i].coordinates[m]);
	while(double_string != NULL)
	{
	m++;
	if(m>=DIM) break;
	double_string = strtok (NULL, " ");
	p_array[i].coordinates[m] = atof(double_string);
	//printf("%.6f ", p_array[i].coordinates[m]);
	}


	m=0;



	//printf("\n");
	//printf("size: %d aaa %d\n",size, i);
	i++;

	}

	printf("number of pairs from generator: %d\n",i);
	fclose(myfile);   //close the file prior to exiting the routine


}

void fill_projected_point_array()
{

	int i=0;

	Point origin;

	for(i=0;i<DIM;i++)
	origin.coordinates[i]=0;

	i=0;
	while(i<size)
	{
			   convert_to_uvw(i);


			  p_array[i].dist_projection_to_origine=calculate_distance(origin,p_array[i]);

			 // printf("\n( %.4f  %.4f  %.4f ) =>  ",arr[i].coordinates[0],arr[i].coordinates[1],arr[i].coordinates[2]);
            //  printf("projected: [%d]= ( %.4f  %.4f  %.4f ) dist: %.4f\n",i,arr[i].coordinates_projected[0],arr[i].coordinates_projected[1],arr[i].coordinates_projected[2],arr[i].dist_projection_to_origine);
			  i++;
	}

}


void random_points()
{
	int i=0; double t=0; int j=0; int n=0; int index_=0;
	char str[150]; int k=0; int m=0;
	char double_[50];

	if((myfile = fopen(return_path_, "w"))==NULL) {
		printf("Cannot opeen file.\n");

	}

	n= size*DIM;
	srand((unsigned int)time(NULL));


	AA[0] =(double)rand()/(M_PI *100);



	for (i=0;i<n;i++)
	{
		t=(double)rand()/(M_PI *100 );

		for (j=0;j<i;j++)

		{
			if (t==AA[j])
			{ i--;

			break;
			}
			else
			{
				AA[i]=t;
			}
		}//end for

	}  //end for

	//

	i=0;
	for(j=0;j<size*DIM;j++)
	{
		sprintf (str, "");

		for(i=0;i<DIM;i++)
		{	 sprintf(double_,"%.8f",AA[j+i]);
		      strcat(str,double_);
		       if (i==DIM-1)
		      {
			   strcat(str,"\n");
			  break;
	        	}
	      	else
			strcat(str," ");
		}
         fprintf(myfile, "%s", str);

		for(m=0;m<DIM;m++)
		{
			p_array[k].coordinates[m]=AA[j+m];
		//	 printf("%.6f ",p_array[k].coordinates[m]);
		}


		// printf("\n");

		 // printf("%s \n",str);
		k++;
		j=j+DIM-1;
	}


//	printf("number of distinct doubles: %d\n",j);
//	print_array(p_array);

	fclose ( myfile ) ;
	free(AA);

}

void Histogram_projected(int flag_)
{
int i = 0;int numvalues = size ; int index=0;   int K=6;
 double maxval = 0 ; double minval=0; double freqsize =0;
 double hash_product=0;  double _diff=0;

 int *frequency_coordinate;


 for (i = 0 ; i < numvalues ; i++)
 {
	 if (p_array[i].coordinates_projected[flag_] > maxval)
      maxval = p_array[i].coordinates_projected[flag_] ;
 }



 minval = maxval ;

 for (i = 0 ; i < numvalues ; i++)
 {
     if (p_array[i].coordinates_projected[flag_]< minval)
      minval = p_array[i].coordinates_projected[flag_] ;
 }

 freqsize = maxval - minval ;

 hash_product= numvalues / freqsize;

 frequency_coordinate=  (int *) malloc(sizeof(int)* (numvalues+1)) ;
 for (i = 0 ; i < numvalues ; i++)
 {
      frequency_coordinate[i] = 0 ;
 }

 for (i = 0 ; i < numvalues ; i++)
 {
      _diff = p_array[i].coordinates_projected[flag_] - minval ;
       index= (int)(_diff * hash_product);
// printf("index: %d diff: %.3f hash_product: %.3f\n",index,_diff,hash_product) ;
 frequency_coordinate[index]++ ;
 }

/*

 for (i = 0 ; i < numvalues ; i++)
 {
      printf("%2d\t|   %d \n",i, frequency_coordinate[i]) ;
 }
 printf("\n") ;

 */

 if(flag2__==0)
 {
 for(i=0;i<numvalues;i++)
 {   // printf("logf(size) : %f , logf(2): %f bolumu: %f \n", logf(size),  logf(2), ( logf(size)/logf(2)) );
       if(frequency_coordinate[i]>= ( logf(size)*K/logf(2)) )
       {
		   printf("\nDIM=%d . koordinatta cizgisel birikme var. i: %d freq. coor: %d\n\n",flag_,i,frequency_coordinate[i]);
         break;
       }

 }


 }




 free(frequency_coordinate);


}
void Histogram(int flag_)
{
 int i = 0;int numvalues = size ; int index=0;   int K=6;
 double maxval = 0 ; double minval=0; double freqsize =0;
 double hash_product=0;  double _diff=0;

 int *frequency_coordinate;


 for (i = 0 ; i < numvalues ; i++)
 {
	 if (p_array[i].coordinates[flag_] > maxval)
      maxval = p_array[i].coordinates[flag_] ;
 }



 minval = maxval ;

 for (i = 0 ; i < numvalues ; i++)
 {
     if (p_array[i].coordinates[flag_]< minval)
      minval = p_array[i].coordinates[flag_] ;
 }

 freqsize = maxval - minval ;

 hash_product= numvalues / freqsize;

 frequency_coordinate=  (int *) malloc(sizeof(int)* (numvalues+1)) ;
 for (i = 0 ; i < numvalues ; i++)
 {
      frequency_coordinate[i] = 0 ;
 }

 for (i = 0 ; i < numvalues ; i++)
 {
      _diff = p_array[i].coordinates[flag_] - minval ;
       index= (int)(_diff * hash_product);
// printf("index: %d diff: %.3f hash_product: %.3f\n",index,_diff,hash_product) ;
 frequency_coordinate[index]++ ;
 }

/*

 for (i = 0 ; i < numvalues ; i++)
 {
      printf("%2d\t|   %d \n",i, frequency_coordinate[i]) ;
 }
 printf("\n") ;

 */

 if(flag2__==0)
 {
 for(i=0;i<numvalues;i++)
 {   // printf("logf(size) : %f , logf(2): %f bolumu: %f \n", logf(size),  logf(2), ( logf(size)/logf(2)) );
       if(frequency_coordinate[i]>= ( logf(size)*K/logf(2)) )
       {

        // printf("\nDIM=%d . koordinatta cizgisel birikme var. i: %d freq. coor: %d\n\n",flag_,i,frequency_coordinate[i]);
		 break;
       }

 }


 }




 free(frequency_coordinate);




}




void merge_sort(Point arr[],int low,int high,int sort_type_flag,int coordinate_no)
{
	int mid;
	if(low<high) {
		mid=(low+high)/2;
		// Divide and Conquer
		merge_sort(arr,low,mid,sort_type_flag, coordinate_no);
		merge_sort(arr,mid+1,high,sort_type_flag,coordinate_no);
		// Combine
		merge(arr,low,mid,high,sort_type_flag,coordinate_no);
	}


}

void merge(Point arr[],int l,int m,int h,int sort_type_flag,int coordinate_no)
{
	int n1; int n2; int i; int j; int k;


	n1=m-l+1;
	n2=h-m;

	for(i=0; i<n1; i++)
	{

		arr1[i]=arr[l+i];
	}

	for(j=0; j<n2; j++)
	{
		arr2[j]=arr[m+j+1];

	}

	switch(sort_type_flag)
	{
	case 1:

		arr1[i].coordinates[coordinate_no]=(double) INT_MAX;  // To mark the end of each temporary array
		arr2[j].coordinates[coordinate_no]=(double)INT_MAX;

		i=0;
		j=0;
		for(k=l; k<=h; k++)
		{ //process of combining two sorted arrays
			if(arr1[i].coordinates[coordinate_no]<=arr2[j].coordinates[coordinate_no])
			{
				arr[k]=arr1[i];

				i++;
			}
			else
			{
				arr[k]=arr2[j];

				j++;
			}
		}

		break;


	case 3:


		arr1[i].coordinates_projected[coordinate_no]=(double) INT_MAX;  // To mark the end of each temporary array
		arr2[j].coordinates_projected[coordinate_no]=(double)INT_MAX;

		i=0;
		j=0;
		for(k=l; k<=h; k++)
		{ //process of combining two sorted arrays
			if(arr1[i].coordinates_projected[coordinate_no]<=arr2[j].coordinates_projected[coordinate_no])
			{
				arr[k]=arr1[i];

				i++;
			}
			else
			{
				arr[k]=arr2[j];

				j++;
			}
		}

		break;




	}


	}


	Result bruteForceClosestPair(Point p_array_[],int size_,int type)
	{

	double r;
	Point p={0};
	Point q={0};
	Result result;
	int n=0;

	double Dist;
	int i=0;int j=0;



	r =(double) INT_MAX;
	//result.min_distance=(double)INT_MAX;

	if(size_>1)
	{
	for( i = 0; i <=size_-1; i++)
	{
		for( j = i+1; j <=size_-1; j++)
		{


			Dist= calculate_distance(p_array_[i],p_array_[j]);

			if(Dist < r && Dist!=0)
			{
				r = Dist;
				p= p_array_[i];
				q = p_array_[j];


			}

		}


	}

	result.p=p;
	result.q=q;
	result.min_distance= r;
	}
	else
	{
	result.p=p;
	result.q=q;
	r =(double) INT_MAX;
	result.min_distance=(double)INT_MAX;

	}





	return result;

	}

	double calculate_distance(Point p,Point q)
	{
		double d=0; double d_=0; int i=0;
		for(i=0;i<DIM;i++)
			d_ = d_+ ( (p.coordinates[i]-q.coordinates[i]) * (p.coordinates[i]-q.coordinates[i]) ) ;

		d=sqrtf(d_);
		return d;
	}

	Result smartForce(Point *A_[],int flag)
	{
	double d_min=(double)INT_MAX; double d; double coordinates_min[DIM];  int n=size ,r,i;
	Point closest1, closest2; Result result; int g=0; double total_coordinates_min=0;



	for(r=1; r<=n-1; r++ )
	{
		total_coordinates_min=0;
		for(g=0;g<DIM;g++)
				 coordinates_min[g]=(double)INT_MAX;

		for(i=0;i<n-r;i++)
		{



			for(g=0;g<DIM;g++)
			{
					d=calculate_distance(A_[g][i+r],A_[g][i]);
					if(d!=0 )
					{
					  if(flag==1)
				{
					if((A_[g][i+r].coordinates[g]-A_[g][i].coordinates[g])<coordinates_min[g])
						coordinates_min[g]=(A_[g][i+r].coordinates[g]-A_[g][i].coordinates[g]);
				}

				if(flag==3 )
				{



					if((A_[g][i+r].coordinates_projected[g]-A_[g][i].coordinates_projected[g])<coordinates_min[g])
						coordinates_min[g]=(A_[g][i+r].coordinates_projected[g]-A_[g][i].coordinates_projected[g]);


				}

				if(d<d_min)
				{
					d_min=d;

					closest1= A_[g][i+r];
					closest2= A_[g][i];


				}
					}
					else
						continue;

			}


		}

		for(g=0;g<DIM;g++)
			total_coordinates_min= total_coordinates_min + coordinates_min[g] * coordinates_min[g];
		total_coordinates_min=sqrtf(total_coordinates_min);

		if( d_min <= total_coordinates_min )
			break;
	}

	result.p=closest1;
	result.q=closest2;
	result.min_distance=d_min;
	return result;

	}

	Result recursiveClosestPair(Point A_[], int _size)
	{

	Result r_;
	Result result;
	Result Left;
	Result Right;
	Result closest;
	int N;
	int ns_sequences={0};
	double xm; int i=0;int j=0; int k=0;int ns=0; double _dist;int u=0;int m=0; int h=0; int dim=0;


	int size_=0;

	Point *Left_coordinates={0};
	Point *Right_coordinates={0};
	Point *S_coordinates={0};

	r_.min_distance=(double)INT_MAX;
	result.min_distance=(double)INT_MAX;
	Left.min_distance=(double)INT_MAX;
	Right.min_distance=(double)INT_MAX;
	closest.min_distance=(double)INT_MAX;
	N= _size;
   // printf("N: %d\n",N);





	if(N<=3)
	{
		result= bruteForceClosestPair(A_,N,1);
		//counting++;
		//if(closest.min_distance==0)
			//printf ("counting : %d\n", counting);

		//print_array_with_size(A_[0],N);

		return result;
	}
	else
	{

				Right_coordinates= (Point *) malloc(sizeof(Point)* N);
				Left_coordinates= (Point *) malloc(sizeof(Point)* N);



		   if(N%2==0)
		   {
		   size_=N/2;
		   array_copy(Left_coordinates, A_,size_,0);
		    array_copy(Right_coordinates, A_,N,size_);
		   }
		   else
		   {
		   size_=(N/2)+1;
		   array_copy(Left_coordinates, A_,size_,0);
		   array_copy(Right_coordinates, A_,N,size_-1);

		   }

			//size_ = ceilf(N/2);


		   /*
         counting++;
		 if(counting==73731)
			 scanf("%d",&a);
	     printf("size_: %d  N= %d counting: %d\n",size_,N,counting);
			*/

		/*
		printf("xL:\n");
		print_array(Left_coordinates[0]);



		 printf("xR:\n");
		 print_array(Right_coordinates[0]);
		 */

		xm=A_[size_-1].coordinates[0];


           for(i=0;i<N;i++)
		     {
				 if(A_[i].coordinates[0]<=xm)
	           		{
	           			Left_coordinates[u]=A_[i];

	           			u++;

	           		}
		           	else
		           	{
		           		Right_coordinates[m]=A_[i];

		           		m++;

		           	}

	           	}



		/*
		printf("yL:\n");
		print_array(Left_coordinates[1]);



		 printf("yR:\n");
		 print_array(Right_coordinates[1]);

		 printf("zL:\n");
		print_array(Left_coordinates[2]);



		 printf("zR:\n");
		 print_array(Right_coordinates[2]);
		 */

		Left=recursiveClosestPair(Left_coordinates,size_);
		Right=recursiveClosestPair(Right_coordinates,size_);
		result= Right;


		if(Left.min_distance<Right.min_distance)
		{
			result= Left;
		}



				//free(S_coordinates[i]);
				free(Right_coordinates);
				free(Left_coordinates);




				S_coordinates= (Point *) malloc(sizeof(Point)* N);
				//Right_coordinates= (Point *) malloc(sizeof(Point)* N);
				//Left_coordinates= (Point *) malloc(sizeof(Point)* N);





		    for(i=0;i<N;i++)
		    {
		    	if(fabsf(xm-A_[i].coordinates[0])<result.min_distance)
			    {
			    	S_coordinates[j]=A_[i];

			    	j++;
		    	}

		    }
			ns_sequences=j;
			j=0;


		closest=result;



			for(i=0;i<ns_sequences-1;i++)
			{
				k=i+1;
				while(k<ns_sequences && fabsf(S_coordinates[k].coordinates[0]-S_coordinates[i].coordinates[0])<closest.min_distance)
				{

					//++;

					_dist=calculate_distance(S_coordinates[k],S_coordinates[i]);
					//if(_dist==0)
					//	printf("counting : %d\n",counting);

					if(_dist<closest.min_distance && _dist!=0 )
					{
						r_.min_distance=_dist;

						r_.p =S_coordinates[k];

						r_.q= S_coordinates[i];



						closest=r_;
					}

					k++;
				}
			}

	//counting++;
		//if(closest.min_distance==0)
			//printf ("counting : %d\n", counting);
		//free(Left_coordinates);
		//free(Right_coordinates);
		//free(S_coordinates);


				free(S_coordinates);
				//free(Right_coordinates);
				//free(Left_coordinates);



		return closest;
	}

	}

