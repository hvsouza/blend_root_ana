#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void write_binary(){
  FILE *fout;

  fout = fopen("test.dat","wb");
  uint16_t val[12] = {0,0,0,0,1,1,1,1,0,0,0,0};

  int Size = sizeof(val)/sizeof(val[0]);
  int aux = 0;
  for(int j=0; j<Size; j++) {
    if(aux < 1) fwrite(&val[j] , 1 , 2, fout);
    else if (aux == 1) aux = -1;
    aux++;
  }
  fclose(fout);
}

void read_binary(){
  FILE *fout;

  fout = fopen("test.dat","rb");
  uint32_t val[6] = {0};
  
  int Size = sizeof(val)/sizeof(val[0]);

  for(int j=0; j<Size; j++) {
    fread((char *) &val[j],2,1,fout);
    cout << val[j] << " ";
  }
  cout << endl;


}

void test_binary(){
  write_binary();
  read_binary();
    

  
}
