function [LRBF_psi1,LRBF_psi2,LRBF_u30,LRBF_W3] = comparison(R,Me)
switch R/Me
    case 0.4
        LRBF_psi1=load('psi1_c20_RM0.4_NN441_h0.012.txt');
end 

switch R/Me
    case 0.4
        LRBF_psi2=load('psi2_c20_RM0.4_NN441_h0.012.txt');
end

switch R/Me
    case 0.4
        LRBF_u30=load('u30_c20_RM0.4_NN441_h0.012.txt');
end   

switch R/Me
    case 0.4
        LRBF_W3=load('w3_c20_RM0.4_NN441_h0.012.txt');
end
  

end
