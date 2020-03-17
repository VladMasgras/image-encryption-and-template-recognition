#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
void xorshift(unsigned int w,unsigned int h,unsigned int **R,unsigned int seed)
{
    unsigned int r,k;
    r=seed;
    *R=malloc(sizeof(unsigned int)*(2*w*h-1));
    for (k=1;k<2*w*h;k++)
    {
        r=r^r<<13;
        r=r^r>>17;
        r=r^r<<5;
        (*R)[k]=r;
    }
}
void permutare(unsigned int **p,unsigned int w,unsigned int h,unsigned int *R)
{
    unsigned int n=w*h,r,aux;
    int k;
    *p=malloc(sizeof(unsigned int)*n);
    for (k=0;k<n;k++)
        (*p)[k]=k;
    for (k=n-1;k>=1;k--)
    {
        r=R[n-k]%(k+1);
        aux=(*p)[r];
        (*p)[r]=(*p)[k];
        (*p)[k]=aux;
    }
}
void liniarizare(char* nume_imagine,unsigned int **v,unsigned int *latime_img,unsigned int *inaltime_img,unsigned char **header,int *padding)
{
    FILE* fin=fopen(nume_imagine,"rb");
    if(fin == NULL)
   	{
   		printf("nu am gasit imaginea");
   		return;
   	}
    unsigned int i,j;
    fseek(fin,18,SEEK_SET);
    fread(&(*latime_img),sizeof(unsigned int),1,fin);
    fread(&(*inaltime_img),sizeof(unsigned int),1,fin);
    *v=malloc(sizeof(unsigned int)*(*latime_img)*(*inaltime_img));
    fseek(fin,0,SEEK_SET);
    if((*latime_img) % 4 != 0)
        *padding = 4 - (3 * (*latime_img)) % 4;
    else
        *padding = 0;
    unsigned int k=0;
    *header = malloc(54 * sizeof(unsigned char));
    fread(*header, sizeof(char), 54, fin);
    fseek(fin,0,SEEK_END);
    for (i=0;i<*inaltime_img;i++)
    {

        fseek(fin,-3*(*latime_img)-(*padding),SEEK_CUR);
        for(j=0;j<(*latime_img);j++)
        {
            fread(&(*v)[k],3,1,fin);
            k++;
        }
        fseek(fin,-3*(*latime_img),SEEK_CUR);
    }
    fclose(fin);
}
void deliniarizare(char* nume_fisier,unsigned int **v,unsigned int latime_img,unsigned int inaltime_img,unsigned char **header,int padding)
{
    FILE* fout=fopen(nume_fisier,"wb");
    int i,j,k;
    fwrite(*header, sizeof(unsigned char), 54, fout);
    for (i=inaltime_img-1;i>=0;i--)
    {
        for (j=0;j<latime_img;j++)
        {
            fwrite(&(*v)[latime_img*i+j],3,1,fout);
        }
        for (k=0;k<padding;k++)
            fputc(0,fout);
    }
    fclose(fout);
}
void criptare(char* imagine_sursa,char* imagine_criptata,char* fsecret)
{
   	FILE* fkey=fopen(fsecret,"r");
   	if(fkey == NULL)
   	{
   		printf("nu am gasit fisierul cu key");
   		return;
   	}
    unsigned int latime_img, inaltime_img,n,i,j,*v;
    unsigned char *header;
    unsigned int key,sv,*R;
    int padding;
    fscanf(fkey,"%u%u",&key,&sv);
    liniarizare(imagine_sursa,&v,&latime_img,&inaltime_img,&header,&padding);
    xorshift(latime_img,inaltime_img,&R,key);
    unsigned int *p;
    permutare(&p,latime_img,inaltime_img,R);
    n=latime_img*inaltime_img;
    unsigned int *pp=malloc(sizeof(unsigned int)*n);
    for (i=0;i<n;i++)
        pp[p[i]]=v[i];
    unsigned int *c=malloc(sizeof(unsigned int)*n);
    c[0]=sv^pp[0]^R[n];
    for (i=1;i<n;i++)
    {
        c[i]=c[i-1]^pp[i]^R[n+i];
    }
    deliniarizare(imagine_criptata,&c,latime_img,inaltime_img,&header,padding);
    free(c);
    free(v);
    free(R);
    free(p);
    free(pp);
}
void invers(unsigned int **p,unsigned int n,unsigned int **sigma)
{
    for (unsigned int i=0;i<n;i++)
    (*sigma)[(*p)[i]]=i;
}
void decriptare(char* imagine_criptata,char* imagine_decriptata,char* fsecret)
{
    FILE* fkey=fopen(fsecret,"r");
   	if(fkey == NULL)
   	{
   		printf("nu am gasit fisierul cu key");
   		return;
   	}
    unsigned int latime_img, inaltime_img,n,i;
    unsigned char *header;
    unsigned int key,sv,*R;
    unsigned int *c;
    int padding;
    liniarizare(imagine_criptata,&c,&latime_img,&inaltime_img,&header,&padding);
    n=latime_img*inaltime_img;
    fscanf(fkey,"%u%u",&key,&sv);
    xorshift(latime_img,inaltime_img,&R,key);
    unsigned int *p;
    permutare(&p,latime_img,inaltime_img,R);
    unsigned int *sigma=malloc(sizeof(unsigned int)*n);
    invers(&p,n,&sigma);
    unsigned int *cp=malloc(sizeof(unsigned int)*n);
    cp[0]=sv^c[0]^R[n];
    for (i=1;i<n;i++)
    {
        cp[i]=c[i-1]^c[i]^R[n+i];
    }
    unsigned int *d=malloc(sizeof(unsigned int)*n);
    for (i=0;i<n;i++)
        {
            d[sigma[i]]=cp[i];
        }
    deliniarizare(imagine_decriptata,&d,latime_img,inaltime_img,&header,padding);
    free(cp);
    free(c);
    free(R);
    free(p);
    free(d);
    free(sigma);
}
void chi(char* imagine)
{
    unsigned int latime_img, inaltime_img,*c,i;
    unsigned char *header;
    int padding;
   	liniarizare(imagine,&c,&latime_img,&inaltime_img,&header,&padding);
   	unsigned int n=latime_img*inaltime_img;
   	double fm=latime_img*inaltime_img/256,chiR=0,chiG=0,chiB=0;
    unsigned int *vR,*vG,*vB;
    vR=malloc(sizeof(unsigned int)*256);
    vG=malloc(sizeof(unsigned int)*256);
    vB=malloc(sizeof(unsigned int)*256);
    for (i=0;i<=255;i++)
    {
        vB[i]=0;
        vG[i]=0;
        vR[i]=0;
    }
    for (i=0;i<n;i++)
    {
        vB[(c[i] >> (8*0)) & 0xff]++;
        vG[(c[i] >> (8*1)) & 0xff]++;
        vR[(c[i] >> (8*2)) & 0xff]++;
    }
    for(i=0;i<=255;i++)
    {
        chiB+=((vB[i]-fm)*(vB[i]-fm))/fm;
        chiG+=((vG[i]-fm)*(vG[i]-fm))/fm;
        chiR+=((vR[i]-fm)*(vR[i]-fm))/fm;
    }
    printf("Rosu %5.2lf \nVerde %5.2lf \nALbastru %5.2lf\n",chiR,chiG,chiB);
    free(c);
}
void grayscale_image(char* nume_fisier_sursa,char* nume_fisier_destinatie)
{
   FILE *fin, *fout;
   unsigned int dim_img, latime_img, inaltime_img;
   unsigned char pRGB[3], header[54], aux;

   printf("nume_fisier_sursa = %s \n",nume_fisier_sursa);

   fin = fopen(nume_fisier_sursa, "rb");
   if(fin == NULL)
   	{
   		printf("nu am gasit imaginea sursa din care citesc");
   		return;
   	}

   fout = fopen(nume_fisier_destinatie, "wb+");

   fseek(fin, 2, SEEK_SET);
   fread(&dim_img, sizeof(unsigned int), 1, fin);
   printf("Dimensiunea imaginii in octeti: %u\n", dim_img);

   fseek(fin, 18, SEEK_SET);
   fread(&latime_img, sizeof(unsigned int), 1, fin);
   fread(&inaltime_img, sizeof(unsigned int), 1, fin);
   printf("Dimensiunea imaginii in pixeli (latime x inaltime): %u x %u\n",latime_img, inaltime_img);

   //copiaza octet cu octet imaginea initiala in cea noua
	fseek(fin,0,SEEK_SET);
	unsigned char c;
	while(fread(&c,1,1,fin)==1)
	{
		fwrite(&c,1,1,fout);
		fflush(fout);
	}
	fclose(fin);

	//calculam padding-ul pentru o linie
	int padding;
    if(latime_img % 4 != 0)
        padding = 4 - (3 * latime_img) % 4;
    else
        padding = 0;

    printf("padding = %d \n",padding);

	fseek(fout, 54, SEEK_SET);
	int i,j;
	for(i = 0; i < inaltime_img; i++)
	{
		for(j = 0; j < latime_img; j++)
		{
			//citesc culorile pixelului
			fread(pRGB, 3, 1, fout);
			//fac conversia in pixel gri
			aux = 0.299*pRGB[2] + 0.587*pRGB[1] + 0.114*pRGB[0];
			pRGB[0] = pRGB[1] = pRGB[2] = aux;
        	fseek(fout, -3, SEEK_CUR);
        	fwrite(pRGB, 3, 1, fout);
        	fflush(fout);
		}
		fseek(fout,padding,SEEK_CUR);
	}
	fclose(fout);
}
int MAX(int a,int b)
{
    if (a>b) return a;
    return b;
}
int MIN(int a,int b)
{
    if (a<b) return a;
    return b;
}
struct detectie
{
    unsigned char r,g,b;
    int x,y;
    double corr;
};
int compare (const void * a, const void * b)
{

    struct detectie va = *(struct detectie *)a;
    struct detectie vb = *(struct detectie *)b;
    if(va.corr<vb.corr)
        return 1;
    if(va.corr>vb.corr)
        return -1;
    return 0;
}
unsigned int** liniarizare_mat(char* nume_imagine,unsigned int *latime_img,unsigned int *inaltime_img,unsigned char **header,unsigned int *padding)
{
    FILE* fin=fopen(nume_imagine,"rb");
    if(fin == NULL)
    {
        printf("nu am gasit imaginea");
        return;
    }
    unsigned int i,j;
    fseek(fin,18,SEEK_SET);
    fread(&(*latime_img),sizeof(unsigned int),1,fin);
    fread(&(*inaltime_img),sizeof(unsigned int),1,fin);
    unsigned int **v =malloc(*inaltime_img * sizeof(unsigned int *));
    //printf("%u %u",*latime_img,*inaltime_img);
    for (i=0; i<*inaltime_img; i++)
        v[i] =malloc((*latime_img) * sizeof(unsigned int));
    fseek(fin,0,SEEK_SET);
    if((*latime_img) % 4 != 0)
        *padding = 4 - (3 * (*latime_img)) % 4;
    else
        *padding = 0;
    *header = malloc(54 * sizeof(unsigned char));
    fread(*header, sizeof(char), 54, fin);
    fseek(fin,0,SEEK_END);
    for (i=0; i<*inaltime_img; i++)
    {

        fseek(fin,-3*(*latime_img)-(*padding),SEEK_CUR);
        for(j=0; j<(*latime_img); j++)
        {
            fread(&v[i][j],3,1,fin);
        }
        fseek(fin,-3*(*latime_img),SEEK_CUR);
    }
    fclose(fin);
    return v;
}
void deliniarizare_mat(char* nume_fisier,unsigned int **v,unsigned int latime_img,unsigned int inaltime_img,unsigned char **header,int padding)
{
    FILE* fout=fopen(nume_fisier,"wb");
    int i,j,k;
    fwrite(*header, sizeof(unsigned char), 54, fout);
    for (i=inaltime_img-1; i>=0; i--)
    {
        for (j=0; j<latime_img; j++)
        {
            fwrite(&v[i][j],3,1,fout);
        }
        for (k=0; k<padding; k++)
            fputc(0,fout);
    }
    fclose(fout);
}
void colorare(char* nume_imagine,int j,int i,int latime_sablon,int inaltime_sablon,unsigned char r,unsigned char g,unsigned char b)
{
    unsigned int latime_imagine,inaltime_imagine;
    int padding_imagine,k,l;
    unsigned char *header;
    unsigned int **a=liniarizare_mat(nume_imagine,&latime_imagine,&inaltime_imagine,&header,&padding_imagine);
    unsigned int prgb;
    prgb=prgb&0;
    prgb=prgb|r;
    prgb=prgb<<8;
    prgb=prgb|g;
    prgb=prgb<<8;
    prgb=prgb|b;
    for (k=i; k<i+inaltime_sablon; k++)
        for (l=j; l<j+latime_sablon; l++)
        {
            a[i][l]=prgb;
            a[i+inaltime_sablon-1][l]=prgb;
            a[k][j]=prgb;
            a[k][j+latime_sablon-1]=prgb;
        }
    deliniarizare_mat(nume_imagine,a,latime_imagine,inaltime_imagine,&header,padding_imagine);
    free(a);
}

void template_matching(char* nume_imagine,char* nume_sablon,int *s,struct detectie **d,unsigned char r,unsigned char g,unsigned char bl,double ps)
{
    unsigned int**a,**b;
    unsigned int latime_img,inaltime_img,padding_imagine,padding_sablon,latime_sablon,inaltime_sablon;
    unsigned char *header_imagine,*header_sablon;
    a=liniarizare_mat(nume_imagine,&latime_img,&inaltime_img,&header_imagine,&padding_imagine);
    b=liniarizare_mat(nume_sablon,&latime_sablon,&inaltime_sablon,&header_sablon,&padding_sablon);
    unsigned int n=inaltime_sablon*latime_sablon;
    int k,l,di=-1,dj,i,j;
    double sm=0,fim;
    for (i=0; i<inaltime_sablon; i++)
        for (j=0; j<latime_sablon; j++)
            sm+=(b[i][j] & 0xff);
    sm=sm/n;
    double sigs=0,sigfi,corr=0;
    double sumsigs=0,sumsigfi;
    for (i=0; i<inaltime_sablon; i++)
        for (j=0; j<latime_sablon; j++)
            sumsigs+=((b[i][j] & 0xff)-sm)*((b[i][j] & 0xff)-sm);
    sigs=sqrt(sumsigs/(n-1));
    for (i=0; i<=inaltime_img-inaltime_sablon; i++)
    {
        di++;
        dj=-1;
        for (j=0; j<=latime_img-latime_sablon; j++)
        {
            dj++;
            sumsigfi=0;
            sigfi=0;
            fim=0;
            corr=0;
            for (k=i; k<i+inaltime_sablon; k++)
                for (l=j; l<j+latime_sablon; l++)
                {
                    fim+=(a[k][l]&0xff);
                }
            fim=fim/n;
            for (k=i; k<i+inaltime_sablon; k++)
                for (l=j; l<j+latime_sablon; l++)
                    sumsigfi+=((a[k][l]&0xff)-fim)*((a[k][l]&0xff)-fim);
            sigfi=sqrt(sumsigfi/(n-1));
            for (k=i; k<i+inaltime_sablon; k++)
                for (l=j; l<j+latime_sablon; l++)
                    corr+=((double)((a[k][l]&0xff)-fim)*((b[k-di][l-dj]&0xff)-sm))/(sigfi*sigs);
            corr=corr/n;
            if (corr>ps)
            {
                (*d)[*s].corr=corr;
                (*d)[*s].x=j;
                (*d)[*s].y=i;
                (*d)[*s].r=r;
                (*d)[*s].g=g;
                (*d)[*s].b=bl;
                (*s)++;
                *d=realloc(*d,((*s)+1)*sizeof(struct detectie));
            }
        }
    }
    for(i=0;i<inaltime_img;i++)
        free(a[i]);
    free(a);
    for(i=0;i<inaltime_sablon;i++)
    free(b[i]);
    free(b);
}
void sortare(struct detectie **d,int s)
{
    qsort(*d,s,sizeof(struct detectie),compare);
}
void eliminare(struct detectie **d,int *s,char*nume_sablon)
{
    int i,j;
    FILE* fin=fopen(nume_sablon,"rb");
    if(fin == NULL)
    {
        printf("nu am gasit sablonul");
        return;
    }
    int latime_sablon,inaltime_sablon,k;
    double si,su,r;
    fseek(fin,18,SEEK_SET);
    fread(&latime_sablon,sizeof(unsigned int),1,fin);
    fread(&inaltime_sablon,sizeof(unsigned int),1,fin);
    fclose(fin);
    for (i=0; i<(*s)-1; i++)
        for (j=i+1; j<(*s); j++)
        {
            si=MAX(0,MIN((*d)[i].x+latime_sablon,(*d)[j].x+latime_sablon)-MAX((*d)[i].x,(*d)[j].x))*MAX(0,MIN((*d)[i].y+inaltime_sablon,(*d)[j].y+inaltime_sablon)-MAX((*d)[i].y,(*d)[j].y));
            su=latime_sablon*inaltime_sablon+latime_sablon*inaltime_sablon-si;
            r=si/su;
            if(r>0.2)
            {
                for (k=j; k<(*s)-1; k++)
                {
                    (*d)[k]=(*d)[k+1];
                }
                (*s)--;
                j--;
            }
        }
}
int main()
{
    char imagine_originala[101],imagine_criptata[101],fisier_secret[101],imagine_decriptata[101];
    printf("Numele fisierului care contine imaginea originala de criptat: ");
    fgets(imagine_originala, 101, stdin);
    imagine_originala[strlen(imagine_originala) - 1] = '\0';
    printf("Numele fisierului care contine imaginea criptata: ");
    fgets(imagine_criptata, 101, stdin);
    imagine_criptata[strlen(imagine_criptata) - 1] = '\0';
    printf("Numele fisierului care contine fiserul cu cheia secreta: ");
    fgets(fisier_secret, 101, stdin);
    fisier_secret[strlen(fisier_secret) - 1] = '\0';
    printf("Numele fisierului care contine imaginea decriptata: ");
    fgets(imagine_decriptata, 101, stdin);
    imagine_decriptata[strlen(imagine_decriptata) - 1] = '\0';
    criptare(imagine_originala,imagine_criptata,fisier_secret);
    decriptare(imagine_criptata,imagine_decriptata,fisier_secret);
    chi(imagine_originala);
    chi(imagine_criptata);

    struct detectie *d;
    d=malloc(sizeof(struct detectie));
    unsigned int latime_sablon,inaltime_sablon,padding,**v;
    unsigned char *header;
    unsigned char r,g,b;
    int s=0,i,n,canalr,canalg,canalb;
    char nume_img_initiala[101],nume_img_initiala_grey[101],nume_sablon[101],nume_sablon_gr[101];
    printf("Cititi imaginea initiala: ");
    fgets(nume_img_initiala, 101, stdin);
    nume_img_initiala[strlen(nume_img_initiala) - 1] = '\0';
    printf("Numele imaginitii initiala in greyscale ");
    fgets(nume_img_initiala_grey, 101, stdin);
    nume_img_initiala_grey[strlen(nume_img_initiala_grey) - 1] = '\0';
    grayscale_image(nume_img_initiala,nume_img_initiala_grey);
    printf("Dati numarul de sabloane: ");
    scanf("%d",&n);
    fgetc(stdin);
    for (i=0;i<n;i++)
        {
            printf("Cititi numele sablonului ");
            fgets(nume_sablon, 101, stdin);
            nume_sablon[strlen(nume_sablon) - 1] = '\0';
            printf("Numele sablonului grayscale ");
            fgets(nume_sablon_gr, 101, stdin);
            nume_sablon_gr[strlen(nume_sablon_gr) - 1] = '\0';
            grayscale_image(nume_sablon,nume_sablon_gr);
            printf("Dati valorile R,G,B despartite printr-un spatiu ");
            scanf("%d %d %d",&canalr,&canalg,&canalb);
            fgetc(stdin);
            r=(unsigned char)canalr;
            g=(unsigned char)canalg;
            b=(unsigned char)canalb;
            template_matching(nume_img_initiala_grey,nume_sablon_gr,&s,&d,r,g,b,0.50);
        }
    sortare(&d,s);
    eliminare(&d,&s,nume_sablon);
    FILE* fin=fopen(nume_sablon,"rb");
    if (fin==NULL)
    {
        printf("Nu s-a gasit sablonul");
        return;
    }
    fseek(fin,18,SEEK_SET);
    fread(&latime_sablon,sizeof(unsigned int),1,fin);
    fread(&inaltime_sablon,sizeof(unsigned int),1,fin);
    fclose(fin);
    for (i=0;i<s;i++)
        colorare(nume_img_initiala,d[i].x,d[i].y,latime_sablon,inaltime_sablon,d[i].r,d[i].g,d[i].b);
    printf("Am terminat de colorat imaginea");
    free(d);
    return 0;
}

