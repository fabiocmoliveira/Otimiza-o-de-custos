function [alpha]=dicotomica1D(f,comprimentointervalofinal,epsilon,intervalo,kmax)

%DADOS DE ENTRADA
% f=@(x)x^2+2*x;
% comprimentointervalofinal=1e-11;
% epsilon=comprimentointervalofinal/10;
% intervalo=[-3 5];
% kmax=100;

k=1;
while k<kmax
    
   
    if intervalo(2)-intervalo(1)<comprimentointervalofinal
        
        alpha=sum(intervalo)/2;
        break
    
    else
    
        lambda=sum(intervalo)/2-epsilon;
        mi=sum(intervalo)/2+epsilon; 
        
        if f(lambda)<f(mi)
        
            intervalo(2)=mi;
        
        else
        
            intervalo(1)=lambda;
            
        end
    
        k=k+1;
    
    end
    
end

if k==kmax
    
    disp('Número máximo de iterações atingido!')
end
% else
%     
%     disp(['O ponto de mínimo é : ' num2str(x)])
% 
% end