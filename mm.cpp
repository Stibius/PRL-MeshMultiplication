/*
 * algorithm: Mesh Multiplication 
 * author: Jan Vybíral, xvybir05
 *
 */

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
 
using namespace std;

constexpr int TAG = 0;
 
/*
Precte soubor se vstupni matici.
inputFile - nazev vstupniho souboru
matrix - sem se ulozi nactena matice
count - sem se ulozi hodnota z prvniho radku souboru
rows - sem se ulozi pocet nactenych radku matice
cols - sem se ulozi pocet nactenych sloupcu matice
navratova hodnota - false, pokud se soubor nepodarilo otevrit nebo obsahuje neplatnou hodnotu 
nebo je pocet nactenych hodnot 0 nebo vsechny radky nemaji stejny pocet hodnot, jinak true
*/
bool readMatrixFile(const string& inputFile, vector<vector<int>>& matrix, int& count, int& rows, int& cols)
{
    string sValue;   //sem se nacita obsah mezi oddelovaci ' ' a '\n'

    fstream fin;     //vstupni stream

    fin.open(inputFile, ios::in);   //otevre vstupni soubor pro cteni  
    //kontrola, zda se soubor podarilo otevrit
    if (!fin.is_open())
    {
        return false;
    }

    char c;            //sem se postupne nacitaji jednotlive znaky ze souboru
    bool first = true; //indikuje, ze jsme na prvnim radku souboru 

    //cti znaky az do konce souboru
    do
    {
        fin.get(c); //nacti jeden znak ze souboru

        //pokud je aktualni znak oddelovac nebo jsme narazili na konec souboru
        if (c == ' ' || c == '\n' || !fin.good())
        {
            //pokud mame neco nacteno
            if (sValue.size() != 0)
            {
                //pokus se to prevest na cele cislo
                try
                {
                    int value = stoi(sValue); //prevod retezce na cele cislo
                    //pokud jsme na prvnim radku souboru
                    if (first)
                    {
                        count = value; //uloz hodnotu z prvniho radku
                        first = false; //u dalsich hodnot uz se predpoklada, ze nejsou na prvnim radku
                    }
                    //pokud nejsme na prvnim radku souboru
                    else
                    {
                        matrix.back().push_back(value); //pridej danou hodnotu na posledni misto posledniho radku vystupu
                    }
                    if (c == '\n') matrix.push_back(vector<int>()); //pokud jsme na konci radku, pridej novy prazdny radek do vystupu 
                    sValue.clear(); //vymazeme retezec pro dalsi nacitani
                }
                //nepodarilo se prevest retezec na cele cislo
                catch (const invalid_argument & ia)
                {
                    return false;
                }
                //hodnota je mimo rozsah datoveho typu int
                catch (const out_of_range & oor)
                {
                    return false;
                }
            }
        }
        //pokud aktualni znak neni oddelovac ani konec souboru
        else
        {
            sValue += c; //nacti znak hodnoty
        }
    } while (fin.good());

    fin.close(); //zavri vstupni soubor

    //zjisti pocet nactenych radku a zkontroluj, jestli neni 0
    rows = matrix.size();
    if (rows == 0)
    {
        return false;
    }

    //na konci muze byt prazdny radek, pokud tam je, odstran ho
    if (matrix.back().size() == 0)
    {
        matrix.erase(matrix.end() - 1);
    }

    //zjisti pocet zbylych radku a zkontroluj, jestli neni 0
    rows = matrix.size();
    if (rows == 0)
    {
        return false;
    }

    //zjisti pocet nactenych sloupcu u prvniho radku a zkontroluj, jestli neni 0
    cols = matrix[0].size();
    if (cols == 0)
    {
        return false;
    }
    //zkontroluj, jestli maji vsechny radky stejny pocet hodnot
    for (int i = 1; i < matrix.size(); ++i)
    {
        if (matrix[i].size() != cols)
        {
            return false;
        }
    }

    return true;
}
 
/*
Precte oba vstupni soubory.
matrix1 - sem se ulozi nactena matice ze souboru "mat1"
matrix2 - sem se ulozi nactena matice ze souboru "mat2"
M = sem se ulozi pocet radku prvni matice
K - sem se ulozi pocet sloupcu druhe matice
N - sem se ulozi pocet sloupcu/radku prvni/druhe matice
navratova hodnota - true, pokud se vse podarilo nacist a data jsou korektni, jinak false
*/
bool readInputs(vector<vector<int>>& matrix1, vector<vector<int>>& matrix2, int& M, int& K, int& N)
{
    int count;  //sem se ulozi hodnota z prvniho radku souboru 
    int N1;     //sem se ulozi pocet sloupcu prvni matice
    int N2;     //sem se ulozi pocet radku druhe matice

    //nacti prvni matici a zkontroluj, jestli pocet jejich radku odpovida hodnote z prvniho radku souboru
    if (!readMatrixFile("mat1", matrix1, count, M, N1) || count != M)
    {
        return false;
    }

    //nacti druhou matici a zkontroluj, jestli pocet jejich sloupcu odpovida hodnote z prvniho radku souboru
    if (!readMatrixFile("mat2", matrix2, count, N2, K) || count != K)
    {
        return false;
    }

    //zkontroluj, jestli je pocet sloupcu prvni matice roven poctu radku druhe matice a lze je tak nasobit
    if (N1 != N2)
    {
        return false;
    }

    N = N1; //uloz hodnotu poctu sloupcu/radku prvni/druhe matice

    return true;
}

int main(int argc, char* argv[])
{
    int numprocs;               //pocet procesoru
    int myid;                   //muj rank
    int a;                      //hodnota, ktera prijde od leveho souseda a preda se pravemu
    int b;                      //hodnota, ktera prijde od horniho souseda a preda se spodnimu
    int c;                      //vysledna hodnota tohoto procesoru
    int M;                      //pocet radku prvni matice
    int N;                      //pocet sloupcu/radku prvni/druhe matice 
    int K;                      //pocet sloupcu druhe matice
    int i;                      //radek aktualniho procesoru
    int j;                      //sloupec aktualniho procesoru

    vector<vector<int>> matrix1;  //sem si nulty procesor nacte prvni matici
    vector<vector<int>> matrix2;  //sem si nulty procesor nacte druhou matici
    vector<int> A(N);             //pokud se jedna o procesor z prvniho sloupce, tak se sem ulozi korespondujici radek z prvni vstupni matice
    vector<int> B(N);             //pokud se jedna o procesor z prvniho radku, tak se sem ulozi korespondujici sloupec ze druhe vstupni matice
    int result;                   //sem nulty procesor prijima vysledky od ostatnich
    MPI_Status stat;              //struct- obsahuje kod- source, tag, error

    //MPI INIT
    MPI_Init(&argc, &argv);                         //inicializace MPI 
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);       //zjistíme, kolik procesů běží 
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);           //zjistíme id svého procesu 

    bool inputOK = true; //indikuje, jestli nedoslo k chybe pri nacitani vstupu
    //nulty procesor nacte a ulozi si obe matice a udaje o jejich rozmerech
    if (myid == 0)
    {
        if (!readInputs(matrix1, matrix2, M, K, N))
        {
            inputOK = false;
        }
    }

    //nulty procesor posle ostatnim informaci o tom, jestli byl vstup v poradku
    MPI_Bcast(&inputOK, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);

    //pokud vstup nebyl v poradku, vsechny procesory se ukonci
    if (inputOK == false)
    {
        MPI_Finalize();

        return -1;
    }

    //nulty procesor posle ostatnim informace o rozmerech matic
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&K, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //vypocet radku a sloupce aktualniho procesoru
    i = myid / K;
    j = myid % K;

    //nulty procesor
    if (myid == 0)
    {
        //ulozi si prvni radek prvni matice
        for (int col = 0; col < N; ++col)
        {
            A[col] = matrix1[0][col];
        }

        //ulozi si prvni sloupec druhe matice
        for (int row = 0; row < N; ++row)
        {
            B[row] = matrix2[row][0];
        }

        //rozesle procesorum v prvnim sloupci (krome sebe) jejich korespondujici radky z prvni matice
        vector<int> A_send(N);
        for (int row = 1; row < M; ++row)
        {
            for (int col = 0; col < N; ++col)
            {
                A_send[col] = matrix1[row][col];
            }
            MPI_Send(A_send.data(), N, MPI_INT, row * K, TAG, MPI_COMM_WORLD);
        }

        //rozesle procesorum v prvnim radku (krome sebe) jejich korespondujici sloupce z druhe matice
        vector<int> B_send(N);
        for (int col = 1; col < K; ++col)
        {
            for (int row = 0; row < N; ++row)
            {
                B_send[row] = matrix2[row][col];
            }
            MPI_Send(B_send.data(), N, MPI_INT, col, TAG, MPI_COMM_WORLD);
        }
    }
    //procesory z prvniho radku krome nulteho procesoru
    else if (i == 0)
    {
        //prijmou od nulteho procesoru svoje korespondujici sloupce z druhe matice
        MPI_Recv(B.data(), N, MPI_INT, 0, TAG, MPI_COMM_WORLD, &stat);
    }
    //procesory z prvniho sloupce krome nulteho procesoru
    else if (j == 0)
    {
        //prijmou od nulteho procesoru svoje korespondujici radky z prvni matice
        MPI_Recv(A.data(), N, MPI_INT, 0, TAG, MPI_COMM_WORLD, &stat);
    }

    //vsechny procesory vynuluji svuj vysledek
    c = 0;

    //provede se tolik iteraci, kolik je sloupcu/radku v prvni/druhe matici
    for (int it = 0; it < N; ++it)
    {
        //procesory v prvnim radku nacitaji hodnoty ze svych vstupnich poli
        if (i == 0)
        {
            b = B[it];
        }
        //procesory v ostatnich radcich prijimaji hodnoty od procesoru o jeden radek nad nimi
        else
        {
            MPI_Recv(&b, 1, MPI_INT, myid - K, TAG, MPI_COMM_WORLD, &stat);
        }

        //procesory v prvnim sloupci nacitaji hodnoty ze svych vstupnich poli
        if (j == 0)
        {
            a = A[it];
        }
        //procesory v ostatnich sloupcich prijimaji hodnoty od procesoru o jeden sloupec vlevo
        else
        {
            MPI_Recv(&a, 1, MPI_INT, myid - 1, TAG, MPI_COMM_WORLD, &stat);
        }

        //postupny vypocet vysledne hodnoty pro tento procesor
        c += a * b;

        //vsechny procesory krome tech v poslednim sloupci preposlou svuj vstup procesorum o jeden sloupec vpravo
        if (j < K - 1) MPI_Send(&a, 1, MPI_INT, myid + 1, TAG, MPI_COMM_WORLD);

        //vsechny procesory krome tech v poslednim radku preposlou svuj vstup procesorum o jeden radek pod nimi
        if (i < M - 1) MPI_Send(&b, 1, MPI_INT, myid + K, TAG, MPI_COMM_WORLD);
    }

    //vsechny procesory krome nulteho poslou svuj vysledek nultemu procesoru
    if (myid != 0)
    {
        MPI_Send(&c, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD);
    }
    //nulty procesor provede vypis cele vystupni matice
    else
    {
        cout << M << ":" << K << "\n";
        for (int row = 0; row < M; ++row)
        {
            for (int col = 0; col < K; ++col)
            {
                //nulty procesor vypise svuj vysledek
                if (row == 0 && col == 0)
                {
                    if (col == K - 1)
                        cout << c << "\n";
                    else
                        cout << c << " ";
                }
                //nulty procesor prijme vysledek jineho procesoru a vypise ho
                else
                {
                    MPI_Recv(&result, 1, MPI_INT, (row * K) + col, TAG, MPI_COMM_WORLD, &stat);

                    if (col == K - 1)
                        cout << result << "\n";
                    else
                        cout << result << " ";
                }
            }
        }
    }

    MPI_Finalize();

    return 0;

}

