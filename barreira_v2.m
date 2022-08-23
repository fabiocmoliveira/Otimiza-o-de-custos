clc
close all
clear all
    
%Função objetivo e suas derivadas
f=@(x)(x(1)-2)^4+(x(1)-2*x(2))^2;
f1=@(x)2*x(1)-4*x(2)+4*(x(1)-2)^3;
f2=@(x)8*x(2)-4*x(1);
fd={f1;f2};
%Restrições e derivadas de B=-1/h
h=@(x)(x(1)^2-x(2));
B=@(x)-1/(x(1)^2-x(2));
B1=@(x)((2*x(1))/(x(2)-x(1)^2)^2);
B2=@(x)(-1/(-x(1)^2+x(2))^2);
Bd={B1;B2};

%Ponto inicial
x0=[0 1]';
    
%Precisão para a convergência
epsilon=1e-5;

%Parâmetros para a busca dicotômica
parametros_dicotomica={epsilon,100,[0 1]};

%Parâmetros para o BFGS
parametros_BFGS={epsilon,100,h};

%Penalização inicial
c=10;

%Incremento para o parêmetro de penalização
alfa=1/10;

%Número máximo de iterações
imax=3;
    
%Contagem das iterações
i=0;
while i<imax

    %Função de penalização
    F=@(x)f(x)+c*B(x);
    %Gradiente da função de penalização
    gradientes=cell(length(x0),1);    
    gradientes{1}=@(x)fd{1}(x)+c*Bd{1}(x);
    GRAD=@(x)gradientes{1}(x);
    for j=2:length(x0)
        
        gradientes{j}=@(x)fd{j}(x)+c*Bd{j}(x);
        GRAD=@(x)[GRAD(x); gradientes{j}(x)];
        
    end
   
%     %Função de penalização
%     F=@(x)f(x)+c*B(x);
%     %Gradiente da função de penalização
%     GRAD=@(x)[f1(x)+c*B1(x);f2(x)+c*B2(x)];
    
    %Minimização irrestrita da função penalizada: método BFGS
    x1=metodo_BFGS(F,GRAD,x0,parametros_BFGS,parametros_dicotomica);

    %Convergência
    if abs(x1-x0)<=epsilon
            
        break
            
    end
    
    %Atualização do parâmetro de penalidade e do ponto
    c=alfa.*c;
    x0=x1;
    i=i+1;
    
end
disp(['numero de iteracoes: ' num2str(i)])
disp('ponto otimo x*:')
disp(num2str(x1))