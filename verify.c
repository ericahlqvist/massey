#include <pari/pari.h>
#include <stdio.h>
#include <time.h>

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#include "misc_functions.h"
#include "find_orbit.h"
#include "ext_and_aut.h"
#include "find_variables.h"
#include "find_basis.h"
#include "test_basis.h"
#include "artin_symbol.h"
#include "test_artin.h"

int
main (int argc, char *argv[])	  
{
    //--------
    clock_t start = clock();
    //--------
    // char prime_str = argv[1];
    // char my_det = argv[2];
    //char swap_str[100];
    

    int my_int;

    int min;
    int sec;
    int msec;
    
    pari_init(4000000000,500000);
    // printf("Initial adress: %ld\n", avma);
    // pari_sp limit = stack_lim(avma, 1);
    
    GEN p = cgeti(DEFAULTPREC);
    GEN s = pol_x(fetch_user_var("s"));
    GEN K, f, p_ClFld_pol;

    p = gp_read_str(argv[1]);
    my_int = atoi(argv[2]);

    // Define K.pol
    f = gsubgs(gsqr(s), my_int);

    // Define base field K
    K = Buchall(f, nf_FORCE, DEFAULTPREC);

    p_ClFld_pol = bnrclassfield(K, p, 0, DEFAULTPREC);

    GEN K_ext_aut   =   my_ext(K, p_ClFld_pol, my_int, s, p, 0);
    GEN LxAbs       =   gel(K_ext_aut, 1);
    GEN LyAbs       =   gel(K_ext_aut, 3);

    int ans1 = bnfcertify(K);
    int ans2 = bnfcertify(LxAbs);
    int ans3 = bnfcertify(LyAbs);

    if (ans1 & ans2 & ans3) {
        printf(ANSI_COLOR_GREEN "GRH OK!\n\n" ANSI_COLOR_RESET);
    }
    else {
        printf(ANSI_COLOR_RED "GRH not ok\n\n" ANSI_COLOR_RESET);
    }

    printf(ANSI_COLOR_GREEN "Done! \n \n" ANSI_COLOR_RESET);

    // Close pari
    pari_close();

    //--------
    // Compute the time the whole program took
    clock_t duration = (clock()-start) / 1000;
    msec = duration%1000000;
    sec = (duration/1000)%60;
    min = duration/60000;

    printf (ANSI_COLOR_YELLOW "Runtime: %d min, %d,%d sec\n\n" ANSI_COLOR_RESET, min, sec, msec);
    
    //-----------
    return 0;
}