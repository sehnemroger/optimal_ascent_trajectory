% a=0;
% tic
% for i=1:1e7
%     a=a+1;
% end
% toc
tic
a=1;
while a<=1e7
    a=a+1;
    if a> 10
        break
    end
end
toc
a
i