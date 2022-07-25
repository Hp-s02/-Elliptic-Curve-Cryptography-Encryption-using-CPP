#include <stdio.h>
#include <conio.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <sys\timeb.h>
#include <cstdlib>
#include <iostream>
#include <time.h>
#include <string.h>

using namespace std;

class ecc
{
public :
	ecc();
	void exchange_key();
	void encrypt_ECC();
	void decrypt_ECC();
    long int prime(long int);
	long double mod(long double, long double);
	void random_point(long int, long double, long double, long int);
	long int inverse_devision(long int, long double);
	void xy_twopointadd(long int, long double, long double, long double, long double, long double);
	void xy_Multiplypoint(long int, long int, long int, long double, long double);

private :
	long double x3,y3;
	long double x_random,y_random;
};
ecc::ecc() {}
struct timeb t;

long int ecc::prime(long int len)
{
	srand(time(NULL));
	long int i;
	long int rand_numb,count,result_devision,check,min,max,prime,cal_prim;
	min=pow(10,len-1);
	max=pow(10,len);
	do
    {
		do
		{
			rand_numb=mod(rand(),max);
		}while(rand_numb<=min || rand_numb<=3);
		for(i=2;i<=sqrt(rand_numb);i++)
		{
			count=rand_numb/i;
			result_devision=floor(count);
			check=result_devision*i;
			if(check==rand_numb)
				break;
			else if(check<rand_numb || rand_numb==2)
				continue;
		}
		cal_prim=rand_numb;
	}while(check==rand_numb);
	prime=cal_prim;
	return prime;
}

long double ecc::mod(long double a,long double b)
{
	if(fmodl(a,b)<0)
		return (fmodl(a,b)+b);
	else
		return fmodl(a,b);
}

void ecc::random_point(long int p, long double a4, long double a6, long int keylength)
{
	long int check_point,zero_one,check_root;
	long double check_x;
	if(keylength<=5)
    {
		do
		{
			x_random=rand()%p+1;
			y_random=rand()%p+1;
			check_point=mod(mod(powl(y_random,2),p)-mod(mod(powl(x_random,2),p)*x_random,p)-mod(a4*x_random,p)-a6,p);
		}while(check_point!=0);
	}
	else if(keylength>5)
	{
		do
		{
			x_random=rand()%p+1;
			check_x=sqrtl(mod(mod(mod(powl(x_random,2),p)*x_random,p)+mod((a4*x_random),p)+a6,p));
			check_root=floor(check_x);
		}while(check_x!=check_root);
		zero_one=rand()%2;
		zero_one;
		y_random = zero_one=1 ?  mod(check_x,p) : mod(-check_x,p);
	}
}

long int ecc::inverse_devision(long int p,long double devision)
{
	long double a1,a2,a3,b1,b2,b3,t1,t2,t3;
	long int q;
	a1=1;
	b1=0;
	a2=0;
	b2=1;
	a3=p;
	b3=mod(devision,p);
	do
    {
		q=floor(a3/b3);
		t1=a1-mod(q*b1,p);
		t2=a2-mod(q*b2,p);
		t3=a3-mod(q*b3,p);

		a1=b1;
		a2=b2;
		a3=b3;
		b1=t1;
		b2=t2;
		b3=t3;
	}while(b3!=1);
	return mod(b2,p);
}

void ecc::xy_twopointadd(long int p,long double a4, long double x1, long double x2, long double y1, long double y2)
{
	long double pembilangx,devisionx,pembilangy,devisiony,inverse_devisionx,inverse_devisiony;
	if(x1==x2 && y1==y2)
	{
		pembilangx=mod(powl(mod(3*mod(powl(x1,2),p)+a4,p),2),p);
		devisionx=mod(powl(mod(2*y1,p),2),p);
		pembilangy=mod(mod(3*mod(pow(x1,2),p)+a4,p),p);
		devisiony=mod(mod(2*y1,p),p);
		if(devisionx!=0)
		{
			inverse_devisionx=inverse_devision(p,devisionx);
			inverse_devisiony=inverse_devision(p,devisiony);
			x3=mod((mod(pembilangx*inverse_devisionx,p))-x1-x2,p);
			y3=mod(mod(pembilangy*inverse_devisiony,p)*mod((x1-x3),p)-y1,p);
		}
		else if(devisionx==0)
		{
			x3=0; // only sign
			y3=0; // only sign
		}
	}
	else if(x1!=x2 || y1!=y2)
	{
		if(x2==0 && y2==0)
		{
			y3=y1;
			x3=x1;
		}
		else if(x1==0 && y1==0)
		{
			y3=y2;
			x3=x2;
		}
		else
		{
			pembilangx=mod(powl(y2-y1,2),p);
			devisionx=mod(powl(x2-x1,2),p);
			pembilangy=mod(y2-y1,p);
			devisiony=mod(x2-x1,p);
			if(devisionx!=0)
			{
				inverse_devisionx=inverse_devision(p,devisionx);
				inverse_devisiony=inverse_devision(p,devisiony);
				x3=mod((mod(pembilangx*inverse_devisionx,p))-x1-x2,p);
				y3=mod(mod(pembilangy*inverse_devisiony,p)*mod((x1-x3),p)-y1,p);
			}
			else if(devisionx==0)
			{
				x3=0; //only sign
				y3=0; //only sign
			}
		}
	}
}

// binary algorithm
void ecc::xy_Multiplypoint(long int p,long int a4, long int k, long double x, long double y)
{
	long int u[100],l,j,x1,y1;
	l=1;
	while(k>0)
    {
		if(k%2==1)
		{
			if(k==1)
			{
				u[l]=1;
				break;
			}
			u[l]=1;
		}
		else
			u[l]=0;
		k=k/2;
		l=l+1;
	}

	// addition with binary
	x3=0;
	y3=0;
	for(j=l;j>0;j--)
    {
		if(u[j]==1)
		{
			x1=x3;
			y1=y3;
			xy_twopointadd(p,a4,x1,x1,y1,y1);
			x1=x3;
			y1=y3;
			xy_twopointadd(p,a4,x1,x,y1,y);
		}
		else if(u[j]==0)
		{
			x1=x3;
			y1=y3;
			xy_twopointadd(p,a4,x1,x1,y1,y1);
		}
	}
}

static int l = 0;
void ecc::exchange_key()
{
	ofstream file_output;
	file_output.open("key.txt");
	ftime(&t);
	long int time_start, time_fthissh;
	float time_gen,second;
	long double a4,a6;
	long int p;
	long int private1,private2;
	long int x_key1,y_key1,x_key2,y_key2;
	long double x_public1,y_public1,x_public2,y_public2;
	char strbit[1000],buffer[1000];
	long int bit,count_input;
	long double det;

	srand(time(NULL));

	cout<<"\n-------------------Exchange key----------------\n";
	cout<<"\nGenerate prime p:\n";
	count_input=0;
	start:
	cout<<"\t=>Enter key length: ";
	strcpy(strbit,gets(buffer));
	if (atoi(strbit)<=0 || atoi(strbit)>9 ||strlen(strbit)!=1){
		cout<<"\n<WARNING> : 'key length must be integer between 1 to 9'\n\n";
		count_input++;
		if(count_input%4==0){
			system("cls");
		}
		goto start;
	}
	bit=atoi(strbit); // bit is keylength
	time_start=time(NULL);
	p=prime(bit);
	cout<<"\n\tp = "<<p<<"\n";

	cout<<"\nGenerate curve : \n";
	do
    {
		a4=rand()%p;
		a6=rand()%p;
	}while(mod((4*mod(mod(powl(a4,2),p)*a4,p))+(27*mod(powl(a6,2),p)),p)==0);
	det=mod((4*mod(mod(powl(a4,2),p)*a4,p))+(27*mod(powl(a6,2),p)),p);
	cout<<"\t=>a4 = "<<long(a4)<<", "<<"a6 = "<<long(a6)<<"\n";
	cout<<"\t=>4*a4^3 + 27*a6^2 (mod "<<p<<") = " <<long(det)<<"\n\n";
	cout<<"Elliptic curve:"<<" "<<"y^2 = x^3 + "<<long(a4)<<"x + "<<long(a6);

	cout<<"\n\nDecide random point (x,y) : \n";
	random_point(p, a4, a6, bit);
	cout<<"\t=>Generated random point = ("<<long(x_random)<<","<<long(y_random)<<")";

	cout<<"\n\nGenerate private key : \n";
	do
    {
		private1=rand()%(p-1)+1;
		private2=rand()%(p-1)+1;
	}while(private1<1 || private2<1);

	cout<<"\tprivate1 = "<<private1<<"\n";
	cout<<"\tprivate2 = "<<private2<<"\n\n";

	cout<<"Calculate public key : \n";
	xy_Multiplypoint(p,a4,private1,x_random,y_random);
	x_public1=x3;
	y_public1=y3;
	xy_Multiplypoint(p,a4,private2,x_random,y_random);
	x_public2=x3;
	y_public2=y3;
	cout<<"\t=>public1 ="<<private1<<"*("<<long(x_random)<<","<<long(y_random)<<")=("<<long(x_public1)<<","<<long(y_public1)<<")\n";
	cout<<"\t=>public2 ="<<private2<<"*("<<long(x_random)<<","<<long(y_random)<<")=("<<long(x_public2)<<","<<long(y_public2)<<")\n\n";

	cout<<"Calculate private key : \n";
	xy_Multiplypoint(p,a4,private1,x_public2,y_public2);
	x_key1=x3;
	y_key1=y3;
	xy_Multiplypoint(p,a4,private2,x_public1,y_public1);
	x_key2=x3;
	y_key2=y3;
	cout<<"\t=>key1 ="<<private1<<"*("<<long(x_public2)<<","<<long(y_public2)<<")=("<<long(x_key1)<<","<<long(y_key1)<<")\n";
	cout<<"\t=>key2 ="<<private2<<"*("<<long(x_public1)<<","<<long(y_public1)<<")=("<<long(x_key2)<<","<<long(y_key2)<<")\n\n";

	file_output<<p<<" "<<long(a4)<<" "<<long(a6)<<" ";
	file_output<<long(x_random)<<" "<<long(y_random)<<" ";
	file_output<<private1<<" "<<private2<<" ";
	file_output<<long(x_public1)<<" "<<long(y_public1)<<" ";
	file_output<<long(x_public2)<<" "<<long(y_public2)<<" ";
	file_output<<x_key1<<" "<<y_key1<<" ";
	file_output<<x_key2<<" "<<y_key2<<" ";
	time_fthissh=time(NULL);
	second=(t.millitm/1000.0);
	time_gen=(time_fthissh-time_start)+second;
	cout<<endl;
	cout<<"The time cost of generating key = "<<time_gen<<" "<<"second"<<endl;
}

void ecc::encrypt_ECC()
{
	ftime(&t);
	long int time_start, time_fthissh;
	float time_gen,second;
	long int x_public1_gen,y_public1_gen;
	long int x_key1_gen;
	long int loop,i,loop_p;
	long int iter;
	long int k,sum_block;
	char fplaintext_in[500],fplaintext_out[500],buffer[500],fplaintext_inTemp[500],fplaintext_outTemp[500];
	long int count_input;

	srand(time(NULL));
	again1:
	cout<<"\n\n-------------------ECC Encryption----------------\n\n";
	cout<<"\nGenerate prime p :\n";

	char datakey[]="key.txt";
	long double *m_key;
	ifstream file_dtkey;

	file_dtkey.open(datakey);
	m_key=new long double[50];
	loop=0;
	while(!file_dtkey.eof())
	{
		file_dtkey>>m_key[loop+1];
		loop++;
	}
	iter=loop-1; //because read by number
	file_dtkey.close();
	if(iter==0)
    {
		cout<<"\n<WARNING> : 'File "<<datakey<<" is NULL, please enter keylength to key exchange!'\n\n";
		exchange_key();
		system("cls");
		goto again1;
	}
	for(i=1;i<=iter;i++)
		cout<<long(m_key[i])<<" ";
	cout<<"\n\n";

	cout<<"User 1's information of Encryption:"<<"\n";
	cout<<"\t=>prime = "<<long(m_key[1])<<"\n";
	cout<<"\t=>a4 = "<<long(m_key[2])<<"\n";
	cout<<"\t=>a6 = "<<long(m_key[3])<<"\n";
	cout<<"\t=>random point(P) = ("<<long(m_key[4])<<","<<long(m_key[5])<<")\n";
	cout<<"\t=>public2 = ("<<long(m_key[10])<<","<<long(m_key[11])<<")\n";
	cout<<"\t=>private1 = "<<long(m_key[6])<<"\n\n";

	count_input=0;
	again:
	strcpy(fplaintext_out,"ciphertext_");

	cout<<"Please enter the file name you want to encrypt = ";
	buffer[0]=char(500);
	strcpy(fplaintext_in,gets(buffer));

	FILE *input;
	long int m_plaintext;
	input = fopen(fplaintext_in, "r");

	if (input== NULL)
    {
		cout<<"\n<WARNING> : 'File "<< fplaintext_in<<" not found, enter name file plaintext with correct!'\n\n";

		count_input++;
		if(count_input%4==0)
			system("cls");
		goto again;
	}

	strcpy(fplaintext_inTemp,strrev(fplaintext_in));
	int pos_slash=0;
	for(i=1;i<=int(strlen(fplaintext_inTemp));i++){
		pos_slash++;
		if(fplaintext_inTemp[i]==92)
			break;
	}
	pos_slash= int(strlen(fplaintext_inTemp)) - pos_slash; // from example key.txt become chipertext_key.txt

	//pengkodisian if pos_slash = 0 or != 0
	if(pos_slash==0)
        // if on current folder
		//kond. 1"<<"\n";
		strcpy(fplaintext_out,strcat(fplaintext_out,strrev(fplaintext_in))); //fplaintext_out = embed fplaintext_in with chipertext_.txt
	else
	{ // if with path dynamic
		//kond. 2"<<"\n";
		strcpy(fplaintext_in,strrev(fplaintext_in));

		// split name path become 2 part :
		//part1
		int loop_take=0;
		for(i=0;i<pos_slash;i++)
        {
			fplaintext_outTemp[loop_take]=fplaintext_in[i];
			loop_take++;
		}
		fplaintext_outTemp[loop_take]= '\0';  // to destroy noise
		strcpy(fplaintext_out,strcat(fplaintext_outTemp,fplaintext_out));

		//part2
		loop_take=0;
		for(i=pos_slash;i< int(strlen(fplaintext_inTemp));i++)
        {
			fplaintext_outTemp[loop_take]=fplaintext_in[i];
			loop_take++;
		}
		fplaintext_outTemp[loop_take]= '\0';  // to destroy noise
		strcpy(fplaintext_out,strcat(fplaintext_out,fplaintext_outTemp));
	}
	cout<<"\nThe resulting file of encryption = ";
	printf(fplaintext_out);
	cout<<"\n";

	time_start=time(NULL);
	ofstream file_output; //file fplaintext_out
	file_output.open(fplaintext_out);

	cout<<"\nprocess encrypt/ user1 Count( plaintext[] ^ x_key1) :"<<"\n";
	cout<<"Please wait.";
	loop_p=0;
	sum_block=1000;
	while((m_plaintext=fgetc(input))!= EOF)
    {
		if(loop_p%sum_block==0)
		{
			do
			{
				k=mod(rand(), m_key[1]); //k is private1_gen have characteristic dynamic
			}while(k<2);

			xy_Multiplypoint(m_key[1],m_key[2],k,m_key[4],m_key[5]); //x_Multiplypoint(prime,a4,private1_gen,xrandom,yrandom)
			x_public1_gen=x3;
			y_public1_gen=y3; // m_key[4] is x_random, m_key[5] is y_random
			file_output<<x_public1_gen<<" "<<y_public1_gen<<" ";

			//k*public key user2
			xy_Multiplypoint(m_key[1],m_key[2],k,m_key[10],m_key[11]); //x_Multiplypoint(prime,a4,private1_gen,xpublic2,ypublic2)
			x_key1_gen=x3;
		}
		m_plaintext=(m_plaintext^x_key1_gen); //XOR
		file_output<<std::hex<<m_plaintext<<" ";
		loop_p++;
		//animation event process encrypt
		if(loop_p%100000==0)
			cout<<".";
        l++;
	}
	cout<<"!";
	iter=loop_p; //because read by character
	file_output<<endl;
	fclose(input);



	cout<<"\n\nsum character plaintext = "<<iter;
	cout<<"\n\n";
	cout<<"Encryption Completed!\n";
	time_fthissh=time(NULL);
	second=(t.millitm/1000.0);
	time_gen=(time_fthissh-time_start)+second;
	cout<<endl;
	cout<<"The time cost of encryption = "<<time_gen<<" "<<"second"<<endl;

	ifstream inFile;
	inFile.open(fplaintext_out);
	string str;

	while(getline(inFile,str))
    {
		cout << str << endl;
	}
}

int main ()
{
	//declaration name object from class
	ecc run;
	char str_choice[1000],buffer[1000];
	buffer[0]=char(1000); // this is space to patch character input start from null
	int _choice;
	cout.setf(ios::fixed,ios::floatfield);   // float field set to fixed
	do
    {
		system("cls");

		printf("===========================================\n");
		printf("\tElliptic Curve Cryptography");
		printf("\n===========================================\n");
		cout<<"\t1. Create exchange key\n";
		cout<<"\t2. Encrypt ECC\n";
		cout<<"\t3. Exit \n";
		cout<<"\nEnter your choice: ";
		strcpy(str_choice,gets(buffer));
		_choice=atoi(str_choice);
		switch(_choice)
		{
            case 1:
            {
                run.exchange_key();
                cout<<endl;
                cout<<"\npress any key to be continue...... \n";
                getch ();
                break;
            }
            case 2:
            {
                run.encrypt_ECC();
                cout<<"\n\npress any key to be continue...... \n";
                getch ();
                break;
            }
		}
	}while (_choice!=3);
}
