#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"

#define M_PI 3.14159265358979323846


/*
Smooths a grayscale image by convolving it with a Gaussian kernel of standard deviation sigma.
Input:
    Image im: the input image
    float sigma: the standard deviation of the Gaussian kernel
Output:
    Image: the smoothed image (im.w, im.h, 1)
*/


//Il primo passo è la riduzione del rumore dell'immagine applicando un filtro Gaussiano.
//Per farlo, utilizzate le funzioni che avete sviluppato nel homework 2, e completate
//la funzione smooth_image().
//Usando un filtro Gaussiano di dimensione 9x9 e deviazione standard σ=1.4
Image smooth_image(const Image& im, float sigma)
{
    // TODO: Your code here
    Image gauss = make_gaussian_filter(sigma);
    Image smooth = convolve_image(im,gauss,false);
    return smooth;
}


/*
Computes the magnitude and direction of the gradient of an image.
Input:
    Image im: the input image
Output:
    pair<Image,Image>: the magnitude and direction of the gradient of the image
                       with magnitude in [0,1] and direction in [-pi,pi]
*/

/*
Il secondo passo consiste nell'identificare l'intensità e la direzione dei bordi nell'immagine. 
Questo si può ottenere calcolando il gradiente dell'immagine. Lo step 1 è necessario perché
i gradienti sono molto sensibili al rumore. Per calcolare il gradiente, completate 
la funzione compute_gradient(). Potete utilizzare la funzione sobel che avete sviluppato 
nel homework 2, ricordando di normalizzare l'ampiezza del gradiente.
*/
pair<Image,Image> compute_gradient(const Image& im)
{
    // TODO: Your code here
    pair<Image,Image> sobel = sobel_image(im);
    feature_normalize(sobel.first);
    return sobel;
}


/*
Performs non-maximum suppression on an image.
Input:
    Image mag: the magnitude of the gradient of the image [0,1]
    Image dir: the direction of the gradient of the image [-pi,pi]
Output:
    Image: the image after non-maximum suppression
*/
/*
L'ampiezza del gradiente ha messo in evidenza i bordi, ma questi sono troppo spessi 
e sfocati. Questo step ha lo scopo di ridurre la larghezza dei bordi e di renderli 
più definiti. Fondamentalmente, questo è fatto preservando i pixel che sono localmente
massimi lungo la direzione del bordo ed eliminando tutti gli altri pixel.
Nel nostro caso, invece di utilizzare direttamente la direzione del gradiente, 
la arrotondiamo al multiplo di 45° più vicino. Poi consideriamo i valori dei due pixel
vicini lungo questa direzione. In pratica, invece di ottenere i valori p ed r per 
interpolazione, stiamo usando l'approccio nearest neighbor.

Più precisamente, per ogni pixel:

Si arrotonda la direzione del gradiente al più vicino multiplo di PI/4.
Se il pixel è maggiore dei suoi vicini lungo questa direzione, allora viene mantenuto, 
altrimenti viene eliminato.
*/
float min_mag(float a, float b, float c)
  {
  return min({a,b,c});
  }

Image non_maximum_suppression(const Image& mag, const Image& dir)
{
    Image nms(mag.w, mag.h, 1);
    float neighbor1, neighbor2;

    // Iterate through the image and perform non-maximum suppression
    for (int y = 0; y < mag.h; y++) {
        for (int x = 0; x < mag.w; x++) {
            
            // TODO: Your code here

            // Get the direction of the gradient at the current pixel
            float curr_dir_f = (float)(dir(x,y,0)/M_PI)*180;
            int curr_dir = round(curr_dir_f);
            if (curr_dir < 0) curr_dir += 180;
            int near_dir;

            // Round the direction to the nearest multiple of PI/4
            int dist = curr_dir % 45;
            if (dist > 22) near_dir = curr_dir - dist + 45;
            else near_dir = curr_dir - dist;
            

            // Get the magnitude of the gradient of the two neighbors along that direction
            // (Hint: use clamped_pixel to avoid going out of bounds)
            
            if (near_dir == 0 || near_dir == 180) {
                neighbor1 = mag.clamped_pixel(x-1,y,0);
                neighbor2 = mag.clamped_pixel(x+1,y,0);
            }
            else if (near_dir == 90) {
                neighbor1 = mag.clamped_pixel(x,y+1,0);
                neighbor2 = mag.clamped_pixel(x,y-1,0);
            }
            else if (near_dir == 135) {
                neighbor1 = mag.clamped_pixel(x+1,y-1,0);
                neighbor2 = mag.clamped_pixel(x-1,y+1,0);
            }
            else if (near_dir == 45) {
                neighbor1 = mag.clamped_pixel(x-1,y-1,0);
                neighbor2 = mag.clamped_pixel(x+1,y+1,0);
            }

            // If the magnitude of the gradient of the current pixel is greater than that of both neighbors,
            // then it is a local maximum
            if (mag(x,y,0) >= neighbor1 && mag(x,y,0) >= neighbor2) nms.set_pixel(x,y,0,mag(x,y,0));
            else nms.set_pixel(x,y,0,0);
        }
    }

    return nms;
}



/*
    Applies double thresholding to an image.
    Input:
        Image im: the input image
        float lowThreshold: the low threshold value
        float highThreshold: the high threshold value
        float strongVal: the value to use for strong edges
        float weakVal: the value to use for weak edges
    Output:
        Image: the thresholded image
*/

/*
Questo step ha l'obiettivo di identificare 3 tipi di pixel: forti, deboli e non rilevanti:
I pixel forti hanno un valore di intensità maggiore di una soglia alta (high threshold) e sono 
considerati come bordi.
I pixel non rilevanti hanno un valore di intensità minore di una soglia bassa (low threshold)
e sono considerati non rilevanti per il bordo.
I pixel deboli hanno un valore di intensità compreso tra le soglie bassa e alta e sono considerati 
come bordi solo se sono connessi a pixel forti, altrimenti sono considerati non rilevanti. 
Questa operazione verrà eseguita nel prossimo step tramite il tracciamento dei bordi per isteresi.
Per implementare questo step, completate la funzione double_thresholding(). Il risultato di questo
step è un'immagine con 3 valori di intensità: 0, 0.25 e 1.0. I pixel con valore 0 sono quelli non rilevanti,
quelli con valore 1.0 sono quelli forti e quelli con valore 0.25 sono quelli deboli:
*/
Image double_thresholding(const Image& im, float lowThreshold, float highThreshold, float strongVal, float weakVal)
{
    Image res(im.w, im.h, im.c);

    // TODO: Your code here
    for (int c=0; c<im.c; c++) {
        for (int y = 0; y < im.h; y++) {
            for (int x = 0; x < im.w; x++) {
                
                float curr_pix = im(x,y,c);
                if (curr_pix > highThreshold) res.set_pixel(x,y,c,strongVal);
                else if(curr_pix < lowThreshold) res.set_pixel(x,y,c,0);
                else res.set_pixel(x,y,c,weakVal);
            }
        }
    }

    return res;
}


/*
    Applies hysteresis thresholding to an image.
    Input:
        Image im: the input image
        float weak: the value of the weak edges
        float strong: the value of the strong edges
    Output:
        Image: the image after hysteresis thresholding, with only strong edges
*/

/*
Lo step successivo consiste nel trasformare i pixel deboli in pixel forti,
se e solo se sono connessi a pixel forti. In caso contrario, questi pixel vengono impostati a 0.
*/
Image edge_tracking(const Image& im, float weak, float strong)
{
    Image res(im.w, im.h, im.c);

    for (int c=0; c<im.c; c++) {
        for (int y=0; y < im.h; ++y) {
            for (int x=0; x < im.w; ++x) {
                // TODO: Your code here
                if (im(x,y,c) == strong) { res.set_pixel(x,y,c,strong); continue; }
                if (im(x,y,c) == 0) { res.set_pixel(x,y,c,0); continue; }
                
                int set = 0;
                for (int i=-1; i<2; i++) {
                    for (int j=-1; j<2; j++) {
                        
                        if(im.clamped_pixel(x+i,y+j,c) == strong) {
                            res.set_pixel(x,y,c,strong); 
                            set = 1; 
                            break;
                        }   
                    }
                }
                if (set == 0) res.set_pixel(x,y,c,0);
                // Hint: use clamped_pixel when checking the neighbors to avoid going out of bounds
            }
        }
    }
    return res;

}
