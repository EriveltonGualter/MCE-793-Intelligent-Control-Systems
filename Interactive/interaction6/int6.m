% Erivelton Gualter
% 26/03/2018

close all

hold on;

rec = [5 5 7 7; 4 6 6 4];
obstacle =  rectangle(...
        'Position',[5, 4, 2, 2],... % 'Curvature',,...
        'LineWidth',2,...
        'FaceColor',[0.2, 0.7, 0.2],...
        'EdgeColor',0.5*[0.2, 0.7, 0.2]); %Darker version of color

plot([1, 8],[1,7],'ok','LineWidth',2);
    
axis([0 12 0 10]);
drawnow


q0 = [1 1];
qf = [8 7];

qdata = q0;
qplot = q0;

plot(qplot(:,1), qplot(:,2),'-o','LineWidth',2,...
                       'MarkerEdgeColor','b',...
                       'MarkerFaceColor','r',...
                       'MarkerSize',4)
%%
k = 0;
while (k<Inf)
    qk = [12*rand(1) 10*rand(1)];
                   
    [mIdx,mD] = knnsearch(qdata,qk);
    qtemp = [qdata; qk];
    
    idx_qplot = find(qplot(:,1) == qtemp(mIdx,1));
    
    parent = qtemp(mIdx,:);
    xl = linspace(parent(1,1),qk(1,1),100);
    yl = linspace(parent(1,2),qk(1,2),100)
    inl = inpolygon(xl, yl, rec(1,:),rec(2,:))

    if sum(inl) > 1
        idx_l = find(inl,1)-1;
        qk = [xl(idx_l) yl(idx_l)]
    end
    qdata = [qdata; qk];
    
    qplot = [qplot(1:idx_qplot,:); qk; qplot(idx_qplot:end,:)];
    plot(qplot(:,1), qplot(:,2),'-b');
   
    if norm(qk-qf) < 0.01
        break
    end
    pause(1/30)
    k = k + 1
end
               

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); % Full screen
% plot(q0t)