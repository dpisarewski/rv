/*********************************************************
  numbers.h
  --------------------------------------------------------
  generated at Wed Jan 14 14:44:08 2015
  by snns2c ( Bernward Kett 1995 ) 
*********************************************************/

extern "C" int numbers(double *in, double *out, int init);

static struct {
  int NoOfInput;    /* Number of Input Units  */
  int NoOfOutput;   /* Number of Output Units */
  int(* propFunc)(double *, double*, int);
} numbersREC = {160,10,numbers};
