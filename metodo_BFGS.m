function x1=metodo_BFGS(f,g,x0,parametros_BFGS,parametros_dicotomica)
%Entradas

%Parâmetros para a busca dicotômica
% ci=1e-5;
% intmaxdico=100;
% intervalodico=[0 10];
ci=parametros_dicotomica{1};
intmaxdico=parametros_dicotomica{2};
intervalodico=parametros_dicotomica{3};

%Parâmetros para o BFGS
% epsilon=1e-5;
% intmax=100;
epsilon=parametros_BFGS{1};
intmax=parametros_BFGS{2};
restricao=parametros_BFGS{3};

%Início do processo iterativo
int=0;
D=eye(length(x0));%Matriz para a primeira estimativa deve ser positiva definida
while int<intmax

    %Direção de descida
    d=-D*g(x0);d=d/norm(d);

    %Minimização para o passo ótimo 
    falpha=@(alpha)f(x0+alpha*d);
    alphamin=dicotomica1D(falpha,ci,ci/100,intervalodico,intmaxdico);
    
    %Próximo ponto
    x1=x0+alphamin*d;
    %Verifica se o próximo ponto satizfaz as restrições. 
    %Se não, o valor de alphamin é reduzido para 90% até que se encontre um ponto dentro da regiao factível.
    jj=0;
    while restricao(x1)>0
        alphamin=.9*alphamin;
        x1=x0+alphamin*d;
        jj=jj+1;
    end
%     disp(jj)
    %Convergência
    if abs(x1-x0)<epsilon
        
%         disp('Convergência')
        break
        
    end
    
    %Atualização da estimativa da Hessiana
    s=x1-x0;
    y=g(x1)-g(x0);
    if s'*y<=0
        
        D=eye(length(x0));
        disp('reset')
    
    else
        
        r=s/(s'*y)-D*y/(y'*D*y);
        D=D+s*s'/(s'*y)-D*(y*y')*D/(y'*D*y)+r*r'*(y'*D*y);
    
    end
    
    %Número de iterações
    int=int+1;
    x0=x1;
    if int==intmax
        
        disp('Número máximo de iterações atingido')
        
        break
        
    end    
    
end