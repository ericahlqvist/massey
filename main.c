
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
    
    clock_t start = clock();
    

    int p_int, my_int, i;

    int min;
    int sec;
    int msec;
    
    pari_init(8000000000,500000);
    // printf("Initial adress: %ld\n", avma);
    // pari_sp limit = stack_lim(avma, 1);
    
    GEN p = cgeti(DEFAULTPREC);
    GEN s = pol_x(fetch_user_var("s"));
    // GEN x = pol_x(fetch_user_var("x"));
    GEN K, f, Kcyc, p_ClFld_pol, D_prime_vect, D;
    
    p = gp_read_str(argv[1]);
    p_int = atoi(argv[1]);
    my_int = atoi(argv[2]);

    D = stoi(-my_int);
    D_prime_vect = gel(factor(D), 1);
    
    // Define K.pol
    // f = quadpoly0(gneg(D), -1);
    // f = gsubstpol(f, x, s);
    // output(f);
    f = gsubgs(gsqr(s), my_int);
    printf("\n");

    // Define base field K
    K = Buchall(f, nf_FORCE, DEFAULTPREC);
    Kcyc = bnf_get_cyc(K);
    p_ClFld_pol = bnrclassfield(K, p, 0, DEFAULTPREC);

    pari_printf("p_int: %d\n\nmy_int: %d\n\nK_cyc: %Ps\n\nK_basis: %Ps\n\n", p_int, my_int, Kcyc, nf_get_zk(bnf_get_nf(K)));

    GEN K_ext_aut   =   my_ext(K, p_ClFld_pol, my_int, s, p, D_prime_vect, 0);
    GEN LxAbs       =   gel(K_ext_aut, 1);
    GEN LxRel       =   gel(K_ext_aut, 2);
    GEN LyAbs       =   gel(K_ext_aut, 3);
    GEN LyRel       =   gel(K_ext_aut, 4);
    GEN sigma_x     =   gel(K_ext_aut, 5);
    GEN sigma_y     =   gel(K_ext_aut, 6);

    GEN Lx_cyc = bnf_get_cyc(LxAbs);
    GEN Ly_cyc = bnf_get_cyc(LyAbs);

    char file_name[100];
    int Dmod8 = -my_int%8;
    int Dmod16 = -my_int%16;
    int mod;

    if (Dmod8 == 3 || Dmod8 == 7) {
        mod = 8;
    }
    else if (Dmod16 == 4 || Dmod16 == 8) {
        mod = 16;
    }
    else {
        printf(ANSI_COLOR_RED "Wrong discriminant\n\n" ANSI_COLOR_RESET);
        pari_close();
        exit(0);
    }

    sprintf(file_name, "output/%d_%dmod%d.txt", p_int, Dmod8, mod);
    printf("%s", file_name);
    printf("\n");

    FILE *fptr;
    fptr = fopen(file_name, "a");

    pari_fprintf(fptr, "{\"p\": \"%d\", \"D\": \"%d\", \"K-cyc\": \"%Ps\", \"Lx-cyc\": \"%Ps\", \"Ly-cyc\": \"%Ps\", \"Z-rk\": \"-\", \"ZM\": \"-\"},\n", p_int, my_int, Kcyc, Lx_cyc, Ly_cyc);

    fclose(fptr);

    // pari_close();
    // exit(0);

    GEN J_vect = my_find_p_gens(K, p);

    GEN T_x = rnfisnorminit(K, rnf_get_pol(LxRel), 1);
    GEN T_y = rnfisnorminit(K, rnf_get_pol(LyRel), 1);
    
    printf("initializers T_x and T_y found\n\n");

    GEN x_basis = my_find_basis_2(LxAbs, LxRel, K, sigma_x, p, J_vect, T_x);
    GEN y_basis = my_find_basis_2(LyAbs, LyRel, K, sigma_y, p, J_vect, T_y);

    // GEN x_basis = my_find_basis(LxAbs, LxRel, K, sigma_x, J_vect, p);
    // GEN y_basis = my_find_basis(LyAbs, LyRel, K, sigma_y, J_vect, p);

    // Testing basis
    printf("Running basis test ...\n\n");
    my_test_div_a_J_sm1I (LxAbs, LxRel, sigma_x, x_basis);
    my_test_div_b_pI (LxAbs, x_basis, p);
    my_test_sm1_b_pNx_a (LxAbs, LxRel, sigma_x, x_basis, p);

    my_test_div_a_J_sm1I (LyAbs, LyRel, sigma_y, y_basis);
    my_test_div_b_pI (LyAbs, y_basis, p);
    my_test_sm1_b_pNx_a (LyAbs, LyRel, sigma_y, y_basis, p);
    // printf("\nTested if basis in ker(psi)\n\n");

    GEN Ix_vect = my_find_I_from_basis(x_basis);
    // GEN ax_vect = my_find_a_from_basis(x_basis);

    GEN Iy_vect = my_find_I_from_basis(y_basis);
    // GEN ay_vect = my_find_a_from_basis(y_basis);
    
    printf("\n");
    pari_printf("J_vect: %Ps\n\n", J_vect);
    // pari_printf("Ix_vect: %Ps\n\n", Ix_vect);
    // pari_printf("ax_vect: %Ps\n\n", ax_vect);
    // pari_printf("Iy_vect: %Ps\n\n", Ix_vect);
    // pari_printf("ay_vect: %Ps\n\n", ax_vect);
    
    printf("Finding I_prime_vect\n\n");
    GEN Ix_prime_vect = my_find_I_prime_2 (LxAbs, LxRel, K, Ix_vect, sigma_x, p, T_x);
    GEN Iy_prime_vect = my_find_I_prime_2 (LyAbs, LyRel, K, Iy_vect, sigma_y, p, T_y);

    GEN column_1 = zerovec(2);
    GEN column_2 = zerovec(2);

    for (i = 1; i < 3; i++) {
        
        if (p_int == 3)
        {
            
            gel(column_1, i) = stoi(smodis(my_Artin_symbol(LyAbs, LyRel, K, idealmul(K, gel(J_vect, i), rnfidealnormrel(LxRel, gel(Ix_prime_vect, i))), p_int, sigma_y), p_int));

            gel(column_2, i) = stoi(smodis(my_Artin_symbol(LxAbs, LxRel, K, idealmul(K, gel(J_vect, i), rnfidealnormrel(LyRel, gel(Iy_prime_vect, i))), p_int, sigma_x), p_int));
            
        }
        else
        {
            
            gel(column_1, i) = stoi(smodis(my_Artin_symbol(LyAbs, LyRel, K, rnfidealnormrel(LxRel, gel(Ix_prime_vect, i)), p_int, sigma_y), p_int));

            gel(column_2, i) = stoi(smodis(my_Artin_symbol(LxAbs, LxRel, K, rnfidealnormrel(LyRel, gel(Iy_prime_vect, i)), p_int, sigma_x), p_int));
            
        }
        
    }
    
    GEN Zassenhaus_matrix = mkvec2(zerovec(2),zerovec(2));

    
    gmael2(Zassenhaus_matrix, 1,1) = gel(column_1, 1);
    gmael2(Zassenhaus_matrix, 2,1) = gel(column_1, 2);

    gmael2(Zassenhaus_matrix, 1,2) = gel(column_2, 1);
    gmael2(Zassenhaus_matrix, 2,2) = gel(column_2, 2);
    
    printf(ANSI_COLOR_MAGENTA "\n Disc: %d\n\n" ANSI_COLOR_RESET, my_int);

    pari_printf(ANSI_COLOR_MAGENTA "\nPol: %Ps\n\n" ANSI_COLOR_RESET, f);
    
    printf(ANSI_COLOR_YELLOW "Zassenhaus matrix:  \n\n" ANSI_COLOR_RESET);
    pari_printf(ANSI_COLOR_CYAN "%Ps\n\n" ANSI_COLOR_RESET, gel(Zassenhaus_matrix, 1));
    pari_printf(ANSI_COLOR_CYAN "%Ps\n" ANSI_COLOR_RESET, gel(Zassenhaus_matrix, 2));
    printf("\n\n");

    int Z_det = smodis(gsub(gmul(gmael2(Zassenhaus_matrix, 1,1), gmael2(Zassenhaus_matrix, 2,2)), gmul(gmael2(Zassenhaus_matrix, 1,2), gmael2(Zassenhaus_matrix, 2,1))), p_int);

    
    
    FILE    *textfile;
    char    *text;
    long    numbytes;
     
    textfile = fopen(file_name, "r");
     
    fseek(textfile, -26, SEEK_END);
    numbytes = ftell(textfile);
    fseek(textfile, 0L, SEEK_SET);  
 
    text = (char*)calloc(numbytes, sizeof(char));   
 
    fread(text, sizeof(char), numbytes, textfile);
    fclose(textfile);



    fptr = fopen(file_name, "w");
    fprintf(fptr, "%s", text);

    if (my_SQ_MAT_equal0(Zassenhaus_matrix))
    {
        printf(ANSI_COLOR_GREEN "Rank 0  ==> Infinite class field tower\n\n" ANSI_COLOR_RESET);
        pari_fprintf(fptr, " \"Z-rk\": \"0\", \"ZM\": \"%Ps\"},\n", Zassenhaus_matrix);
    }
    else if (Z_det == 0)
    {
        printf(ANSI_COLOR_YELLOW "Rank 1  ==> ZT (3,5), (5,7) or infinite class field tower\n\n" ANSI_COLOR_RESET);
        pari_fprintf(fptr, " \"Z-rk\": \"1\", \"ZM\": \"%Ps\"},\n", Zassenhaus_matrix);
    }
    else {
        printf(ANSI_COLOR_YELLOW "Rank 2  ==> ZT (3,3)\n\n" ANSI_COLOR_RESET);
        pari_fprintf(fptr, " \"Z-rk\": \"2\", \"ZM\": \"%Ps\"},\n", Zassenhaus_matrix);
    }
    
    fclose(fptr);
    
    // output(bnf_get_fu(LxAbs));
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
