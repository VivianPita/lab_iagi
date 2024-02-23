#include "image.h"

// HW0 #1
// const Image& im: input image
// int x,y: pixel coordinates
// int ch: channel of interest
// returns the 0-based location of the pixel value in the data array
int pixel_address(const Image& im, int x, int y, int ch)
  {
  // TODO: calculate and return the index
  
  //NOT_IMPLEMENTED();
  
  return im.w * y + x +(ch*im.w * im.h);
  }

// HW0 #1
// const Image& im: input image
// int x,y,ch: pixel coordinates and channel of interest
// returns the value of the clamped pixel at channel ch
float get_clamped_pixel(const Image& im, int x, int y, int ch)
  {
  // TODO: clamp the coordinates and return the correct pixel value
  
  //NOT_IMPLEMENTED();
  int test_x=x;
  int test_y=y;
  int test_ch=ch;

  if(y<0) test_y=0;
  else if(y>im.h) test_y=-1;

  if(x<0) test_x=0;
  else if(y>im.w) test_x=-1;

  if(ch<0) test_ch=0;
  else if(ch>im.c) test_ch=-1;
  
  return (float) im.data[pixel_address(im,test_x,test_y,test_ch)];
  }


// HW0 #1
// Image& im: input image
// int x,y,ch: pixel coordinates and channel of interest
void set_pixel(Image& im, int x, int y, int c, float value)
  {
  // TODO: Only set the pixel to the value if it's inside the image
  
  //NOT_IMPLEMENTED();
  if(x<0||y<0||c<0||x>=im.w||y>=im.h||c>=im.c) return;
  im(x,y,c)=value;
  }



// HW0 #2
// Copies an image
// Image& to: destination image
// const Image& from: source image
void copy_image(Image& to, const Image& from)
  {
  // allocating data for the new image
  to.data=(float*)calloc(from.w*from.h*from.c,sizeof(float));
  to.c=from.c;
  
  // TODO: populate the remaining fields in 'to' and copy the data
  // You might want to check how 'memcpy' function works
  
  //NOT_IMPLEMENTED();
  to.w=from.w;
  to.h=from.h;
  memcpy(to.data,from.data,to.h*to.w*to.c*(sizeof(float)));
  
  }
