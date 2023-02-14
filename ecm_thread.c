#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <stdlib.h>
#include <pthread.h>
#define NUM_THREADS     90
#define N 1024

typedef struct
{
	mpz_t x;
	mpz_t y;
}
point;

typedef struct thread_data{
	int  thread_id;
	mpz_srcptr B;
	mpz_srcptr x;
	mpz_srcptr y;
	mpz_srcptr A;
	mpz_srcptr n;
}thread_data;

pthread_attr_t attr;
mpz_t n,B;
gmp_randstate_t state;
int num_curves;

int gcd ( int a, int b );
int ec_add(mpz_ptr x3, mpz_ptr y3, mpz_srcptr x1, 
		   mpz_srcptr y1, mpz_srcptr x2, mpz_srcptr y2, 
		   mpz_srcptr A, mpz_srcptr n);
void *fact(void *threadarg);

int gcd ( int a, int b ) 
/*
http://www.math.wustl.edu/~victor/mfmm/compaa/gcd.c
*/
{
	int c;
	while ( a != 0 ) {
		c = a; a = b%a;  b = c;
	}
	return b;
}


int ec_add(mpz_ptr x3, mpz_ptr y3, mpz_srcptr x1, 
		   mpz_srcptr y1, mpz_srcptr x2, mpz_srcptr y2, 
		   mpz_srcptr A, mpz_srcptr n){
	/*
	 (x1,y1) + (x2,y2) 
	 for curve y^2 = x^3 + Ax - B
	 */
	char str[180];
	mpz_t m, m1, m2, b,two, three,negone;
	int success=0, success1;
	mpz_init(two);
	mpz_init(three);
	mpz_init(negone);
	mpz_init(m);
	mpz_init(m1);
	mpz_init(m2);
	mpz_init(b);
	mpz_set_si(two,2);
	mpz_set_si(three,3);
	mpz_set_si(negone,-1);

	if (mpz_cmp(x1,x2) == 0){
		if (mpz_cmp(y1,y2) != 0){
			mpz_clear(m);
		    mpz_clear(m1);
			mpz_clear(m2);
		    mpz_clear(b);
		    mpz_clear(two);
		    mpz_clear(three);
			mpz_clear(negone);
			return 0;
		}
		else{
			success = mpz_invert(m1,two,n);
			success1 = mpz_invert(m2,y1,n);
			mpz_mul(m, m1,m2);
			mpz_mod(m, m ,n); 
			mpz_mul(m1,x1,x1);
			mpz_mod(m1, m1 ,n);
			mpz_mul(m1,three,m1);
			mpz_add(m1,m1,A);
			mpz_mod(m1,m1,n); 
			
			mpz_mul(m,m,m1);
			mpz_mod(m, m ,n);

			if (success == 0 || success1 ==0){
				mpz_gcd(m, y1, n);
				printf("\nFactor: %s\n", mpz_get_str(str,10,m));
				
				mpz_clear(m1);
				mpz_clear(m2);
				mpz_clear(b);
				mpz_clear(two);
				mpz_clear(three);
				mpz_clear(negone);
				
				
				if (mpz_cmp(m,n) == 0){
					mpz_clear(m);
					return -1;
				}
				mpz_clear(m);
				return 0;
			}
		}      
	}
	else{
		mpz_sub(m,x2,x1);
		mpz_mod(m,m,n);
		success = mpz_invert(m,m,n);
		mpz_sub(m1,y2,y1);
		mpz_mul(m,m,m1);
		mpz_mod(m,m,n);

		if (success == 0){
			mpz_sub(m,x2,x1);
			mpz_gcd(m, m, n);	
			printf("\nFactor: %s\n", mpz_get_str(str,10,m));
			
		    mpz_clear(m1);
			mpz_clear(m2);
		    mpz_clear(b);
		    mpz_clear(two);
		    mpz_clear(three);
			mpz_clear(negone);
			if (mpz_cmp(m,n) == 0){
				mpz_clear(m);
				return -1;
			}
			mpz_clear(m);
			return 0;
		}
	}

	mpz_mul(b,m,x1);		   
	mpz_sub(b,y1,b);
	mpz_mod(b,b,n);	   
	mpz_mul(x3,m,m);
	mpz_sub(x3,x3,x1);
	mpz_sub(x3,x3,x2);
	mpz_mod(x3,x3,n);		   
    
	mpz_mul(y3,m,x3);	
	mpz_add(y3,y3,b);
	mpz_mul(y3,y3,negone);	

	mpz_mod(y3,y3,n);	
	mpz_clear(m);
	mpz_clear(m1);
	mpz_clear(m2);
	mpz_clear(b);
	mpz_clear(two);
	mpz_clear(three);
	mpz_clear(negone);
	return 1;					   
}
		
void *fact(void *threadarg){
	// Compute B!(p)
	size_t mystacksize;
    
    char n_string[180];
    point P[33000];
    mpz_t curr_vals[45000];
    mpz_t b_val, curr, g, two,tval;
	mpz_t A,  x,y;
	point pnt;
	pthread_attr_getstacksize (&attr, &mystacksize);
    
	mpz_init(A);
	
	mpz_init(x);
	mpz_init(y);
	
	
	
	mpz_urandomm(A,state, n);
	//mpz_urandomm(x,state, n);
	//mpz_urandomm(y,state, n);
	mpz_set_si(x,0);
	mpz_set_si(y,1);

	mpz_init(b_val);
	mpz_set_si(b_val,2);
	mpz_init(curr);
	mpz_set_si(curr,2);
	mpz_init(g);
	mpz_init(tval);
	mpz_init(two);
	mpz_set_si(two,2);
	
	int m, temp;
    int cursize = 0, psize = 0, i, j;
	
	mpz_init(pnt.x);
	mpz_init(pnt.y);
	mpz_set_si(pnt.x,0);
	mpz_set_si(pnt.y,0);
	
	mpz_init(curr_vals[0]);
    mpz_set(curr_vals[0],curr);
	
	cursize++;
	
	
	m = ec_add(pnt.x, pnt.y, x, y, x, y, A, n);

	
	if (m == 0){
		//free(P);
		mpz_clear(b_val);
		mpz_clear(g);
		mpz_clear(curr);
		mpz_clear(two);
		mpz_clear(pnt.x);
		mpz_clear(pnt.y);
		mpz_clear(x);
		mpz_clear(y);
		pthread_exit(0);
	}
    if (m == -1){
		//free(P);
		mpz_clear(b_val);
		mpz_clear(g);
		mpz_clear(curr);
		mpz_clear(two);
		mpz_clear(pnt.x);
		mpz_clear(pnt.y);
		mpz_clear(x);
		mpz_clear(y);
		pthread_exit(0);
	}
	mpz_set(P[psize].x,pnt.x);  //2P
	mpz_set(P[psize].y,pnt.y);
	mpz_clear(pnt.x);
	mpz_clear(pnt.y);
	psize++;
	for (j=3; mpz_cmp_si(B,j) >= 0; j++){
       
		mpz_set_si(g,j);
		mpz_mul(tval,b_val,g);
		mpz_gcd(g,b_val,g);
		mpz_div(b_val, tval, g);
		
		while(mpz_cmp(curr,b_val) < 0){
			mpz_mul(g,curr,two);
			if (mpz_cmp(g,b_val) <= 0){
				mpz_init(pnt.x);
				mpz_init(pnt.y);
				mpz_set_si(pnt.x,0);
				mpz_set_si(pnt.y,0);
								
				m = ec_add(pnt.x,pnt.y,P[psize-1].x,P[psize-1].y,P[psize-1].x,P[psize-1].y,A,n);
				
				if (m == 0){
					mpz_clear(b_val);
				    mpz_clear(g);
					mpz_clear(curr);
					mpz_clear(two);
					mpz_clear(pnt.x);
				    mpz_clear(pnt.y);
					mpz_clear(x);
					mpz_clear(y);
					pthread_exit(0);
				}
				if (m == -1){
					mpz_clear(b_val);
				    mpz_clear(g);
					mpz_clear(curr);
					mpz_clear(two);
					mpz_clear(pnt.x);
				    mpz_clear(pnt.y);
					mpz_clear(x);
					mpz_clear(y);
					pthread_exit(0);
				}

				mpz_set(P[psize].x,pnt.x);
				mpz_set(P[psize].y,pnt.y);
				mpz_clear(pnt.x);
				mpz_clear(pnt.y);
				psize++;

				mpz_mul(curr, curr, two);
				mpz_init(curr_vals[cursize]);
				mpz_set(curr_vals[cursize],curr);
				cursize++;

			}
			else{
				for (i=cursize-1; i >= 0; i--){

					mpz_add(g,curr,curr_vals[i]);
					if (mpz_cmp(g,b_val) <= 0){
						temp = i;

						mpz_init(pnt.x);
						mpz_init(pnt.y);
						mpz_set_si(pnt.x,0);
						mpz_set_si(pnt.y,0);
						m= ec_add(pnt.x,pnt.y,P[psize-1].x,P[psize-1].y,P[temp].x,P[temp].y,A,n);
						
					
						if (m == 0){
						
							mpz_clear(b_val);
							mpz_clear(g);
							mpz_clear(curr);
							mpz_clear(two);
							mpz_clear(pnt.x);
							mpz_clear(pnt.y);
							mpz_clear(x);
							mpz_clear(y);
							pthread_exit(0);
						}
						if (m == -1){

							mpz_clear(b_val);
							mpz_clear(g);
							mpz_clear(curr);
							mpz_clear(two);
							mpz_clear(pnt.x);
							mpz_clear(pnt.y);
							mpz_clear(x);
							mpz_clear(y);
							pthread_exit(0);
						}
						
						mpz_set(P[psize].x,pnt.x);
						mpz_set(P[psize].y,pnt.y);
						mpz_clear(pnt.x);
						mpz_clear(pnt.y);
						psize++;
						mpz_add(curr,curr,curr_vals[temp]);
						mpz_init(curr_vals[cursize]);
						mpz_set(curr_vals[cursize],curr);
						cursize++;	
					}
				}
			}
		}
	}
	mpz_clear(b_val);
	mpz_clear(g);
	mpz_clear(curr);
	mpz_clear(two);
	mpz_clear(pnt.x);
	mpz_clear(pnt.y);
	mpz_clear(x);
	mpz_clear(y);
	
	
	pthread_exit(0);
}

int main(int argc, char** argv){
	pthread_t threads[NUM_THREADS];
	thread_data thread_data_array[NUM_THREADS];
	char n_string[180];
	int t;
	size_t stacksize;
	
	gmp_randinit_default (state); 
	pthread_attr_init(&attr);
	pthread_attr_getstacksize (&attr, &stacksize);
	stacksize = sizeof(double)*N*N;
	pthread_attr_setstacksize (&attr, stacksize);
	printf("Creating threads with stack size = %li bytes\n",stacksize);
	
	mpz_init(n);
	mpz_init(B);
	
	num_curves = NUM_THREADS;
	
	if(argc == 4){
		num_curves = atoi(argv[3]);
	}
	if (argc < 3){
		mpz_set_si(B,9000);
	}
	else{
		mpz_set_str(B,argv[2], 10);
	}		
	if (argc < 2) {
		mpz_set_si(n, 455839);
		printf("\nUsage: ecm_thread <number> <Bound> <number_curves>\n\n");
		printf("The default Bound is 9000, and default number_curves is 90 \n");
		printf("These are maximum values.  Feel free to specify lower values.");
		printf("Selecting an N for this run.\n\n");
		printf("N is %s\n",mpz_get_str(n_string,10,n));
	}
	else {
		mpz_set_str(n,argv[1], 10);
		
	}
		printf("N is %s\n",mpz_get_str(n_string,10,n));
		printf("%d curves\n",num_curves);
		printf("Bound is %s\n",mpz_get_str(n_string,10,B));
		printf("Use ctrl-c once you have a factor.");
		for(t=0;t<num_curves;t++){
			pthread_create(&threads[t],&attr, fact,(void *) &thread_data_array[t]);//(B,pnt,A,n); 
		}

	pthread_exit(NULL);
}
