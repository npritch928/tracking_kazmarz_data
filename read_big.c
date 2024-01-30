#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int argc, char *argv[])
{
    FILE *ptr,*out;
    char ch;
    int i,stop = 1;
    int count = 0;
    ptr = fopen(argv[1],"r");
    out = fopen(argv[2],"w");
    fprintf(out, "type,value\n rho, ");
    do {
        ch = fgetc(ptr);
        if(ch == '['){
            ++count;
        }else if(count == 3 && ch == ','){ 
            fprintf(out, "\n rho,");
        }else if(count == 4 && ch == ','){
             fprintf(out, "\n iota,");
        }else if(count == 5 && ch == ','){
             fprintf(out, "\n width,");
        }else if(count == 6 && ch == ','){
             fprintf(out, "\n trho,");
        }else if(count == 7 && ch == ']'){
             ch = fgetc(ptr);
	        ++count;
            fprintf(out, "\n iota, ");
            printf("rho\n");
        }else if(count == 4 && ch == ']'){
             ch = fgetc(ptr);
             ++count;
             fprintf(out, "\n width, ");
             printf("iota\n");
        }else if(count == 5 && ch == ']'){
             ch = fgetc(ptr);
             ++count;
             fprintf(out, "\n trho, ");
             printf("width\n");
        }else if(count == 6 && ch == ']'){
             ch = fgetc(ptr);
             ++count;
             printf("trho\ndone\n");
        }else if(count > 2 && stop != 0){
            //if(ch == 'I'){
            //    ptr += 3;
            //}else{
                fprintf(out, "%c",ch); 
            //}
           
        }
        //printf("%d", count);
    }while(ch != EOF);

    fclose(ptr);
    fclose(out);
    return 0;
}
