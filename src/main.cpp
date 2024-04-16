#include <iostream>
#include <stdio.h>

#include "myStats.h"
#include "myRandom.h"

int main()
{
    //int materialSolido[9] = {210, 242, 226, 268, 251, 206, 218, 207};
    
    //showDetails(materialSolido, 9);
    //
    unsigned int len = 365;
    int dias[len], zerados = 0, qtd_commits = 0;
    int dias_desde = 0; // elapsed dias desde 16 de abril de 2024
    bool dias_com_commits[len];

    for(int i = 0; i<len; i++) {
        dias_com_commits[i] = false;
    }

    dias_com_commits[1 + dias_desde] = true;
    dias_com_commits[3 + dias_desde] = true;
    dias_com_commits[4 + dias_desde] = true;
    dias_com_commits[5 + dias_desde] = true;
    dias_com_commits[9 + dias_desde] = true;
    dias_com_commits[14 + dias_desde] = true;
    dias_com_commits[16 + dias_desde] = true;
    dias_com_commits[18 + dias_desde] = true;
    dias_com_commits[19 + dias_desde] = true;
    dias_com_commits[140 + dias_desde] = true;
    dias_com_commits[167 + dias_desde] = true;
    dias_com_commits[175 + dias_desde] = true;
    dias_com_commits[183 + dias_desde] = true;
    dias_com_commits[190 + dias_desde] = true;
    dias_com_commits[191 + dias_desde] = true;
    dias_com_commits[193 + dias_desde] = true;
    dias_com_commits[194 + dias_desde] = true;
    dias_com_commits[195 + dias_desde] = true;
    dias_com_commits[196 + dias_desde] = true;
    dias_com_commits[295 + dias_desde] = true;
    dias_com_commits[296 + dias_desde] = true;

    for (int i=0; i<len; i++) {
        dias[i] = (int) gaussian(2.2, 3.8);
        qtd_commits += dias[i];
        if (dias[i] <= 0) {
            dias[i] = 0;
            zerados++;
        }
    }

    for (int i = 0; i<len; i++) {
        if(dias_com_commits[i]) {
            printf("Dia %i já tem commit\n", i);
            continue;
        }
    }

    for (int i=1; i<=len; i++) {
        if(i >=7 && i <=360)
            continue;

        printf("Dia %i, %i commits realizados\n", i, dias[i-1]);
    }

    showHistogramToScale(dias, len);
    showDetails(dias, len);   
    printf("qtd com zero: %i\ntotal de commits: %i\n", zerados, qtd_commits);
    

    return 0;
}
/*int collection[65536] ;//= {0.7, 0.76, 0.89, 1.3, 0.9, 0.9, 1.1, 1.17, 1.17, 1.2, 1.21, 1.3, 1.3, 1.08, 1.32, 1.32, 1.45, 1.5, 1.6, 1.68, 1.76, 1.82, 1.87, 1.9};
    //char a = 43;
    for(int i=0; i<*(&collection +1)-collection; i++)
    {
        collection[i] = gaussian(146, 4);
    }

    showHistogramToScale(collection, *(&collection +1)-collection);
    showDetails(collection, *(&collection+1)-collection);

    //std::cout << '\\'*true << '\n';
    //putchar(92);*/
    
//#include <cstdlib>
//#include <fstream>
//#include <iostream>
 
//int main()
//{
    //std::system("touch teste/1.txt"); // executes the UNIX command "ls -l >test.txt"
    ////std::cout << std::ifstream("test.txt").rdbuf();
//}
