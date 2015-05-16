/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  main function to launch ellipse center calculation
 *
 *       Compiler:  gcc
 *       vicrucann@gmail.com
 *
 * =====================================================================================
 */
#include "pgm_io.h"
#include "centers.h"

#include <ios>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

image_double average_image(image_double img) {
	int w = img->xsize; 
	int h = img->ysize;
	image_double img_avg = new_image_double_copy(img);
	for (int v = 1; v < h-1; v++) {
		for (int u = 1; u < w-1; u++) {
			double pix = img->data[u-1+(v-1)*w] + img->data[u+(v-1)*w] + img->data[u+1+(v-1)*w] + 
				img->data[u-1+v*w] + img->data[u+v*w] + img->data[u+1+v*w] + 
				img->data[u-1+(v+1)*w] + img->data[u+(v+1)*w] + img->data[u+1+(v+1)*w] ;
			img_avg->data[u+v*w] = pix/9; 
		}
	}
	return img_avg;
}

template <typename T>
int initial_tache(image_double I, vector<T>& h, T& rayon, bool color, T x, T y) {
	int COL_IMA = I->xsize; 
	int LIG_IMA = I->ysize;
	int j = x;
	int i = y;
	int d = 2*rayon;
	if (2*d+1 > LIG_IMA) 
		d=(LIG_IMA-1)/2;
	if(2*d+1 > COL_IMA) 
		d=(COL_IMA-1)/2;
	if(i<d)
		i=d+1;
	if(i>LIG_IMA-1-d) 
		i=LIG_IMA-2-d;
	if(j<d)
		j=d+1;
	if(j>COL_IMA-1-d) 
		j=COL_IMA-2-d;
	int val_haut=0;
	int val_bas=255;
	for (int k = -d; k <= d; k++) {
		for (int l = -d;  l <= d; l++) {
			T lum = I->data[j+l+(i+k)*COL_IMA];
			if (lum > val_haut)
				val_haut = lum;
			else if (lum<val_bas)
				val_bas=lum;
		} }
	T seuil = 0;
	if (!color)
		seuil = val_bas + (val_haut - val_bas)/3 * 2;
	else if (color)
		seuil = val_bas + (val_haut - val_bas)/3;
	
	matrix<T> tab = matrix<T>::zeros(2*d+1, 2*d+1);
	int label = 1;
	
	for (int k = -d+1; k <= d; k++){
		for(int l = -d+1; l <= d-1; l++){
			T lum= I->data[j+l+(i+k)*COL_IMA];
			if (lum < seuil) {
				int imin=l;
				int imax=l+1;          
				while( (I->data[j+imax+(i+k)*COL_IMA] <= seuil) && (imax <= d-1) ){
					imax++; }
				int vallab=0;
				for(int m = imin; m <= imax; m++){
					if (tab(k+d-1,m+d) != 0){
						vallab = tab(k+d-1,m+d);
					}
				}
				if (vallab == 0){
					vallab=label; 
					label++;
				}
				for(int m = imin; m <= imax; m++){
					tab(k+d,m+d)=vallab;
				}
				l=imax;
			}
		}
	}
	matrix<T> bary = matrix<T>::zeros(label, 4);
	for(int k = -d; k <= d; k++){
		for(int l = -d; l <= d; l++){
			if(tab(k+d,l+d)!=0){
				T lum= I->data[j+l+(i+k)*COL_IMA];
				bary(tab(k+d,l+d),0)+=(255-lum)*(j+l);
				bary(tab(k+d,l+d),1)+=(255-lum)*(i+k);
				bary(tab(k+d,l+d),2)+=(255-lum);
				bary(tab(k+d,l+d),3)++; 
			}
		}
	}
	int distmin=100;
	int labelmin=0;
	for(int k = 1; k < label; k++){
		T dist = std::sqrt( (bary(k,0)/bary(k,2)-x) * (bary(k,0)/bary(k,2)-x)+
			(bary(k,1)/bary(k,2)-y) * (bary(k,1)/bary(k,2)-y));
		if(dist < distmin && bary(k,3) > 25 ){ /* 25 =  surface min*/
			distmin=dist;
			labelmin=k;
		}      
	}
	if(labelmin == 0) {
		printf("pb tache trop petite (<=25 pixels)\n");
		return 1;
	}
	x=bary(labelmin,0)/bary(labelmin,2);
	y=bary(labelmin,1)/bary(labelmin,2);
	T sx2 = 0, sy2 = 0, sxy = 0, ss = 0;
	for(int k = -d; k <= d; k++){
		for(int l = -d; l <= d; l++){
			if(tab(k+d,l+d) == labelmin){
				sx2+=(j+l-x)*(j+l-x);
				sy2+=(i+k-y)*(i+k-y);
				sxy+=(i+k-y)*(j+l-x);
				ss++;
			}
		}
	}
	T lambda1 =  ((sx2+sy2)/ss + std::sqrt(( (sx2+sy2)*(sx2+sy2)+4*(sxy*sxy-sx2*sy2)))/ss)/2.0; 
	T lambda2 =  ((sx2*sy2-sxy*sxy)/(ss*ss))/lambda1;
	rayon=std::sqrt(lambda1)*2;
	h[0] = 1.0/rayon; 	                        /* lambda1 		*/
	h[1] = (std::sqrt(lambda1/lambda2))/rayon; 	/* lambda2		*/
	h[2] = std::atan2(sx2/ss-lambda1, -sxy/ss);    	/* alpha  	 	*/
	h[3] = x;        	/* tu     		*/
	h[4] = y;        	/* tv      		*/
	h[5] = 0.25;       	/* rayon cercle 1   	*/
	h[6] = -2.0;      	/* pente 	      	*/
	h[7] = 0.25;       	/* rayon cercle 2   	*/
	h[8] = val_haut; 	/* val_haut   		*/
	h[9] = val_bas; 	/* val_bas   		*/
	h[10] = 1.0; 	/* position step  	*/
	return 0;
}

template <typename T>
vector<T> trgtDataCalc(image_double img_avg, T cx, T cy, T delta) {
	int xbegin = cx+0.5-delta;
	int xend = cx+0.5+delta;
	int ybegin = cy+0.5-delta;
	int yend = cy+0.5+delta;
	int nerr = (yend-ybegin+1)*(xend-xbegin+1);
	vector<T> trgData = vector<T>::zeros(nerr);
	int wi = img_avg->xsize;
	int he = img_avg->ysize;
	int idx = 0;
	for (int v = ybegin; v <= yend; v++) {
		for (int u = xbegin; u <= xend; u++) {
			if (u >= 1 && u <= wi-2 && v >= 1 && v <= he-2)
				trgData[idx] = img_avg->data[u+v*img_avg->xsize]; 
			else
				trgData[idx] = 255; // for 'tache noire'
			idx++;
		}
	}
	return trgData;
}

template <typename T>
T centerLMA(image_double sub_img, bool clr, T& centerX, T& centerY)
{
	image_double img_avg = average_image(sub_img);
	int w = sub_img->xsize; 
	int h = sub_img->ysize;
	T cx = w/2, cy = h/2, radi = 0.4*w;
	vector<T> P(11);
	initial_tache(sub_img, P, radi, clr, cx, cy);
	vector<T> trgData = trgtDataCalc<T>(img_avg, P[3], P[4], radi*2);
	LMTacheC<T> ellipseLMA(img_avg, P[3], P[4], radi*2, clr, w, h);
	T rmse = ellipseLMA.minimize(P, trgData, 0.001);
	free_image_double(img_avg);
    //T lambda1 = P[0]; T lambda2 = P[1]; T theta = P[2];
	centerX = P[3];
	centerY = P[4];
	return rmse;
}

int main(int argc, char ** argv)
{
	bool clr = false; // deals with black circles on white background	
	image_double img = read_pgm_image_double(argv[1]);
	double cx=0, cy=0, l1 = 0, l2 = 0, th = 0;
	centerLMA<double>(img, clr, l1, l2, th, cx, cy);	
free_image_double(img);
	printf("%f %f %f %f %f", l1, l2, th, cx, cy);

	return 0; 	
}
