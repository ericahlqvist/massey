/*
Test basis
*/
void my_test_div_Na_pJ (GEN LxRel, GEN K, GEN my_basis, GEN p) {
    GEN Na, div_Na, pJ, sum;
    int i;
    for (i = 1; i < 3; i++)
    {
        Na = rnfeltnorm(LxRel, gmael2(my_basis, i, 2));
        div_Na = idealhnf(K, Na);
        pJ = idealpow(K, gmael2(my_basis, i, 3), p);
        sum = idealmul(K, div_Na, pJ);
        if (my_SQ_MAT_equal(sum, idealhnf(K, gen_1)))
        {
            printf(ANSI_COLOR_GREEN "Test passed\n\n" ANSI_COLOR_RESET);
        }
        else {
            printf(ANSI_COLOR_RED "my_test_div_Na_pJ failed\n\n" ANSI_COLOR_RESET);
            pari_close();
            exit(0);
        }
    }
    
}


void my_test_div_a_J_sm1I (GEN LxAbs, GEN LxRel, GEN sigma_x, GEN my_basis)
{
    int l = glength(my_basis);
    int i;
    for (i=1; i<l+1; ++i)
    {
        GEN a_elem = gel(gel(my_basis, i),2);
        GEN J_ideal = gel(gel(my_basis, i),3);
        GEN I_ideal = gel(gel(my_basis, i),4);
        GEN div_a_elem = idealhnf(LxAbs, a_elem);

        int equal = my_SQ_MAT_equal(idealmul(LxAbs, div_a_elem, idealmul(LxAbs, rnfidealup0(LxRel, J_ideal, 1), my_1MS_ideal(LxAbs, sigma_x, I_ideal))), idealhnf(LxAbs, gen_1));

        // Print for debugging
        // pari_printf("%Ps\n\n", idealmul(LxAbs, div_a_elem, idealmul(LxAbs, rnfidealup0(LxRel, J_ideal, 1), my_1MS_ideal(LxAbs, sigma_x, I_ideal))));
        // pari_printf("%Ps\n\n",idealhnf(LxAbs, gen_1));


        if (equal)
        {
            printf (ANSI_COLOR_GREEN "Test passed\n" ANSI_COLOR_RESET);
        }
        else
        {
            printf (ANSI_COLOR_RED "ERROR in my_test_div_a_J_sm1I!\n" ANSI_COLOR_RESET);
            pari_printf("%Ps\nNot equal to\n", idealmul(LxAbs, div_a_elem, idealmul(LxAbs, rnfidealup0(LxRel, J_ideal, 1), my_1MS_ideal(LxAbs, sigma_x, I_ideal))));
            pari_printf("%Ps\n\n",idealhnf(LxAbs, gen_1));

            pari_printf("div_a: %Ps\n\nix(J): %Ps\n\n(1-sigma_x)I: %Ps\n\n", div_a_elem, rnfidealup0(LxRel, J_ideal, 1), my_1MS_ideal(LxAbs, sigma_x, I_ideal));

            pari_close();
            exit(0);
        }
    } 
}
/*
Test if div(a) + i_x(J) + (1-\sigma_x)I = 0
*/


void my_test_div_b_pI (GEN LxAbs, GEN my_basis, GEN p) 
{
    int l = glength(my_basis);
    int i;
    for (i=1; i<l+1; ++i)
    {
        GEN b_elem = gel(gel(my_basis, i),1);
        GEN I_ideal = gel(gel(my_basis, i),4);

        GEN div_b_elem = idealhnf(LxAbs, b_elem);
        // pari_printf("%Ps\n\n", idealmul(LxAbs, div_b_elem, idealpow(LxAbs, I_ideal, p)));
        // pari_printf("%Ps\n\n", idealhnf(LxAbs, gen_1));
        if (my_SQ_MAT_equal(idealmul(LxAbs, div_b_elem, idealpow(LxAbs, I_ideal, p)),idealhnf(LxAbs, gen_1)))
        {
            printf (ANSI_COLOR_GREEN "Test passed\n" ANSI_COLOR_RESET);
            // outmat(idealmul(LxAbs, div_b_elem, idealpow(LxAbs, I_ideal, p)));
        }
        else
        {
            printf (ANSI_COLOR_RED "ERROR in my_test_div_b_pI!\n" ANSI_COLOR_RESET);
            pari_printf("%Ps\nNot equal to \n", idealmul(LxAbs, div_b_elem, idealpow(LxAbs, I_ideal, p)));
            pari_printf("%Ps\n\n", idealhnf(LxAbs, gen_1));
            exit(0);
        }
    }
}
/*
Test if (\sigma_x-1)b + (p-N_x)a = 0
*/
void my_test_sm1_b_pNx_a (GEN LxAbs, GEN LxRel, GEN sigma_x, GEN my_basis, GEN p) 
{
    int l = glength(my_basis);
    int i;
    for (i=1; i<l+1; ++i)
    {
        GEN b_elem = gel(gel(my_basis, i),1);
        GEN a_elem = gel(gel(my_basis, i),2);

        GEN elem1 = algtobasis(LxAbs, nfmul(LxAbs, my_SM1_elt(LxAbs, sigma_x, b_elem), my_pmN(LxAbs, LxRel, a_elem, p)));
        
        int equal = my_QV_equal1(elem1);
        
        if (equal)
        {
            printf (ANSI_COLOR_GREEN "Test passed\n" ANSI_COLOR_RESET);
        }
        else
        {
            printf (ANSI_COLOR_RED "ERROR in my_test_sm1_b_pNx_a!\n" ANSI_COLOR_RESET);
            pari_printf("%Ps not equal to \n", elem1);
            pari_printf("%Ps\n\n", algtobasis(LxAbs, gen_1));
            exit(0);
        }
    }
}

int my_test_J_lifts (GEN LxAbs, GEN LxRel, GEN K, GEN J_vect) {

    int l = glength(J_vect);
    GEN J, lift, test_vec;
    int i;
    for (i = 1; i < l+1; i++)
    {
        J = gel(J_vect, i);
        lift = rnfidealup0(LxRel, J, 1);
        test_vec = bnfisprincipal0(LxAbs, lift, 1);
        if (my_QV_equal0(gel(test_vec, 1)))
        {
            printf(ANSI_COLOR_YELLOW "i_x(J[%d]) principal\n\n" ANSI_COLOR_RESET, i);
        }
        else {
            printf(ANSI_COLOR_MAGENTA "i_x(J[%d]) NOT principal\n\n" ANSI_COLOR_RESET, i);
        }
        
    }

    return 0;
}