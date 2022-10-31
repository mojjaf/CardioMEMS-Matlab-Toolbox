function y=Tsallis_entro(x,q)
  
          sum1=sum(x-x.^q);
          sum2=sum1/(q-1);
          y=sum2;
       end