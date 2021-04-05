function x0=step2bCond(f,gradf,grad2f,cond,interval,xi)
% Find zero that makes cond<=0 in interval

E=interval(1); D=interval(2);

n_resolvido=true; % Problema se encontra não resolvido
iter=1; % Iteração atual
pmin=1e-1; % Passo minimo

%%s
% beta_plot=E:.01:D;
% for i=1:length(beta_plot)
%     fplot(i)=f(beta_plot(i));
%     condplot(:,i)=cond(beta_plot(i));
% end
% plot(beta_plot,fplot,'LineWidth',3); hold on
% for i=1:size(condplot,1)
%     plot(beta_plot,condplot(i,:));
% end
%%
while n_resolvido
    x0=step2b(f,gradf,grad2f,xi);
%     %%
%     plot(xi,f(xi),'*');
%     plot(x0,f(x0),'o');
    %%
    if all(cond(x0) <= 0)
        n_resolvido=false; % Sai
    else
        xi_old(iter)=xi; % Salva os chutes iniciais encontrados
%         if abs(xi) > abs(x0) % tende a dar pontos iniciais que puxam para o zero
%             % o proximo zero vai estar depois do próximo minimo ou máximo
%             xi=secT(gradf,grad2f,x0)-pmin;
%         else % zero se encontra a direita do pto inicial
%             xi=secT(gradf,grad2f,x0)+pmin;
%         end
%         if isnan(xi)
%             xi=0;
%         end
%         if xi>D
%             xi=D;
%         elseif xi<E
%                  xi=E;
%         end
        % aparentemente x_util não faz nada, se ele n entra xi n muda.
        % OBSERVAR ISSO 
        xi_util=sum(ismembertol(xi_old,xi,2*pmin));  % verifica valores repetidos
        if xi_util >= 1% Caso o xi ja tenha sido utilizado
            m_inter=0; xi_old_ord=sort(xi_old); xi_old_tam=length(xi_old);
            if xi_old_tam>=10000
                error('Zero não encontrado para as restrições dadas')
            end
            % Escolhe o maior intervalo e a sua localização
            for i=1:xi_old_tam+1
                if i == 1
                    m_inter_cand=xi_old_ord(1)-E;
                elseif i == xi_old_tam+1
                    m_inter_cand=D-xi_old_ord(end);
                else
                    m_inter_cand=xi_old_ord(i)-xi_old_ord(i-1);
                end
                if m_inter_cand>m_inter
                    m_inter=m_inter_cand; %valor do maior intervalo
                    inter=i; % numero do intervalo
                end
            end
            if inter == 1
                xi=E+m_inter/2;
            elseif inter == xi_old_tam % acho que ta errado, deveria ser xi_old_tam+1
                xi=D-m_inter/2;
            else
                xi=m_inter/2+xi_old_ord(inter-1);
            end
        end
    end
    iter=iter+1;
end
end