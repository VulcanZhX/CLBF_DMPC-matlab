function [u] = ustd(Lf,Lg)
    q=sqrt(Lg^2);
    if(q~=0)
       u= -1 * Lg * ((Lf + sqrt(Lf^2 +q^4))/(q^2));
    else
       u = 0;
   end   
end