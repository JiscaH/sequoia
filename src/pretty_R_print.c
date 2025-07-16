#include <R.h>


/* 
* print nicely formatted messages & data from Fortran to R console
*/


/* 
* print timestamp, iteration number, step, no. assigned dams + sires, total likelihood 
* in a table-like format. Print time + iter + step at start of round/step, 
* and no. assigned dams & sires + total likelihood at end of round/step.
*/
void F77_SUB(rprint_status_tbl_header)(void) {
    Rprintf("\n");
    Rprintf("%8s | %2s | %10s | %10s | %5s | %5s | %5s | %10s \n",
            " Time   ", " R", " Step     ", " Progress ", "Dams ", "Sires", " GPs ", "Total LL"); 
    Rprintf("-------- | -- | ---------- | ---------- | ----- | ----- | ----- | ----------\n");
}


void F77_SUB(rprint_status_tbl_entry)(int *time, int *iter, char* step, int *n_parents, double *total_LL) {
    Rprintf("\n");
    Rprintf("%02d:%02d:%02d | %2d | %.10s |            | %5d | %5d | %10.1f \n", 
            time[0], time[1], time[2], *iter, step, n_parents[0], n_parents[1], *total_LL);
}


void F77_SUB(rprint_status_tbl_a)(int *time, int *iter, int *step) {
    char *stepname;
    if (       *step == 1) {
      stepname = "find pairs";
    } else if (*step == 2) {
      stepname = "clustering";
    } else if (*step == 3) {
      stepname = "GP pairs  ";
    } else if (*step == 4) {
      stepname = "merging   ";
    } else if (*step == 5) {
      stepname = "P of sibs ";
    } else if (*step == 6) {
      stepname = "GP Hsibs  ";
    } else if (*step == 7) {
      stepname = "GP Fsibs  ";
    } else if (*step == 8) {
      stepname = "find/check";
    } else if (*step == 90) {
      stepname = "count OH  ";
    } else if (*step == 91) {
      stepname = "est byears";
    } else if (*step == 92) {
      stepname = "Parent LLR";
    } else if (*step == 93) {
      stepname = "Dummy LLR ";
    } else if (*step == 100) {
      stepname = "initial   ";
    } else if (*step == 101) {
      stepname = "ped check ";
    } else if (*step == 102) {
      stepname = "parents   ";
    } else if (*step == 200) {
      stepname = "end       ";
    } else if (*step == 300) {
      stepname = "(all)     ";
    } else {
      stepname = "          ";  
    }
  
    Rprintf("%02d:%02d:%02d | %2d | %.10s | ", time[0], time[1], time[2], *iter, stepname);
}

void F77_SUB(rprint_status_tbl_b)(int *n_parents, int *n_gp, double *total_LL) {
    Rprintf(" | %5d | %5d | %5d | %10.1f \n", n_parents[0], n_parents[1], *n_gp, *total_LL);
}

/* in the for loop of each step, print progress dots at every 10% */
void F77_SUB(rprint_status_tbl_dot)(void) {
    Rprintf(".");
}
 
void F77_SUB(rprint_status_tbl_no_dots)(void) {
    Rprintf("          ");
}

void F77_SUB(rprint_status_tbl_eol)(void) {
    Rprintf(" | \n");
}
 

/* progress bar for GetMaybeRel() etc */
void F77_SUB(rprint_progbar_header)(void) {
    Rprintf(" 0   10  20  30  40  50  60  70  80  90  100%% \n");
    Rprintf(" |   |   |   |   |   |   |   |   |   |   |\n  ");
}

void F77_SUB(rprint_progbar_dot)(void) {
    Rprintf("*");
}

void F77_SUB(rprint_eol)(void) {
    Rprintf("\n");
}