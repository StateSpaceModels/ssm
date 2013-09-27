#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <matheval.h>
     
/* Size of input buffer.  */
#define BUFFER_SIZE 256
     
/* Program is demonstrating use of GNU libmatheval library of procedures
   for evaluating mathematical functions.  */
int
main (int argc, char **argv)
{
char buffer[BUFFER_SIZE];	/* Input buffer.  */
 int length;			/* Length of above buffer. */
 void *f, *f_prim;		/* Evaluators for function and function derivative.  */
 char **names;			/* Function variables names. */
 int count;			/* Number of function variables. */
 double x;			/* Variable x value.  */
 int i;			/* Loop counter. */
     
 /* Read function.  Function has to be over variable x, or result may
    be undetermined.  Size of textual represenatation of function is
    bounded here to 256 characters, in real conditions one should
    probably use GNU readline() instead of fgets() to overcome this
    limit.  */
 printf ("f(x) = ");
 fgets (buffer, BUFFER_SIZE, stdin);
 length = strlen (buffer);
 if (length > 0 && buffer[length - 1] == '\n')
     buffer[length - 1] = '\0';
     
 /* Create evaluator for function.  */
 f = evaluator_create (buffer);
 assert (f);
     
 /* Print variable names appearing in function. */
 evaluator_get_variables (f, &names, &count);
 printf ("  ");
 for (i = 0; i < count; i++)
     printf ("%s ", names[i]);
 printf ("\n");
     
 /* Create evaluator for function derivative and print textual
    representation of derivative.  */
 f_prim = evaluator_derivative_x (f);
 printf ("  f'(x) = %s\n", evaluator_get_string (f_prim));
     
 /* Read variable x value.  */
 printf ("x = ");
 scanf ("%lf", &x);
     
 /* Calculate and print values of function and its derivative for given
    value of x.  */
 printf ("  f(%g) = %g\n", x, evaluator_evaluate_x (f, x));
 printf ("  f'(%g) = %g\n", x, evaluator_evaluate_x (f_prim, x));
     
 /* Destroy evaluators.  */
 evaluator_destroy (f);
 evaluator_destroy (f_prim);
     
 exit (EXIT_SUCCESS);
}
