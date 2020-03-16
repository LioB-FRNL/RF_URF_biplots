function rfps(data,tree,np,prox,pc1,pc2,class,clabels)

% Pseudo-sample representation for random forest models
% Input:
% data - original data matrix
% tree - random forest tree
% np - number of pseudo-samples to be generated in each pseudo-sample
% matrix
% pc1 - index of the first principal component of the data proximity matrix
% to be represented
% pc2 - index of the second principal component of the data proximity matrix
% to be represented
% class - sample class belonging vector
% clabels - cell array containing the class labels
% Output:

%% Setting plotting parameters %%

sym=repmat('sdo*x^+hsdo*x^+h',1,100);
col=repmat('brkymgcbmkbrycbr',1,100);

%% Principal Component Analysis of the double-centered proximity matrix %%

[u,s,v]=svd(doublecentering(1-prox));

disp('PCA model constructed')

%% Generating pseudo-sample matrices %%

psam=zeros(np,size(data,2),size(data,2));

for nvar=1:size(data,2)
    
    ps=zeros(np,size(data,2));
    ps(:,nvar)=linspace(min(data(:,nvar)),max(data(:,nvar)),np);
    psam(:,:,nvar)=ps;
    
end

%% Merging original data and pseudo-sample matrices %%

tree=compact(tree);
tps=zeros(np,size(v,2),size(data,2));
tps_c=zeros(np,size(v,2),size(data,2));

for nvar=1:size(data,2)
    
    mm=[data;squeeze(psam(:,:,nvar))];
    
    %% Estimating pseudo-sample projection scores %%
    
    proxsam=proximity(tree,mm);
    [~,~,~,~,proxsam_c]=doublecentering(1-prox,1-proxsam(size(data,1)+1:end,1:size(data,1)));
    tps(:,:,nvar)=proxsam_c*v;
    tps_c(:,:,nvar)=tps(:,:,nvar)-repmat(mean(tps(:,:,nvar)),size(tps(:,:,nvar),1),1);
    
end

%% Representing original proximity matrix PCA scores plot %%

t=u*s;
evpc1=sum(sum((t(:,pc1)*v(:,pc1)').^2))./(sum(sum((1-prox).^2)));
evpc2=sum(sum((t(:,pc2)*v(:,pc2)').^2))./(sum(sum((1-prox).^2)));

figure

gscatter(t(:,pc1),t(:,pc2),class,col,sym)
axis square
axes=axis;
line([axes(1) axes(2)],[0 0],'LineStyle','--','LineWidth',2,'Color','k')
line([0 0],[axes(3) axes(4)],'LineStyle','--','LineWidth',2,'Color','k')
legend(clabels,'Location','Best')
legend boxoff
xlabel(['PC #',num2str(pc1),' - EV: ',num2str(sprintf('%.2f',evpc1*100))],'FontSize',16,'FontWeight','bold')
ylabel(['PC #',num2str(pc2),' - EV: ',num2str(sprintf('%.2f',evpc2*100))],'FontSize',16,'FontWeight','bold')
box off
set(gca,'FontSize',16,'FontWeight','bold')
set(gcf,'Color','w')

%% Representing pseudo-sample proximity matrix projection scores %%

mt1=abs(max(max(t(:,pc1))));
mt2=abs(max(max(t(:,pc2))));
mtps1=abs(max(max(tps_c(:,pc1,:))));
mtps2=abs(max(max(tps_c(:,pc2,:))));
scf=min(mt1/mtps1,mt2/mtps2);

figure
hold on

for ncl=1:max(class)

    handle(ncl)=scatter(t(class==ncl,pc1),t(class==ncl,pc2),35,class(class==ncl),'filled');
    
end

alpha(handle,0.3);
xlabel(['PC #',num2str(pc1),' - EV: ',num2str(sprintf('%.2f',evpc1*100))],'FontSize',16,'FontWeight','bold')
ylabel(['PC #',num2str(pc2),' - EV: ',num2str(sprintf('%.2f',evpc2*100))],'FontSize',16,'FontWeight','bold')
box off
set(gca,'FontSize',16,'FontWeight','bold')
set(gcf,'Color','w')
cmap_text=colormap(jet(2*size(tps_c,3)));

for nvar=1:size(tps_c,3)
    
    label={[num2str(nvar),'-'];[num2str(nvar),'+']};
    plot(squeeze(tps_c(:,pc1,nvar)).*scf,squeeze(tps_c(:,pc2,nvar)).*scf,'--k')
    text(tps_c([1 end],pc1,nvar).*scf,tps_c([1 end],pc2,nvar).*scf,label,'Color',cmap_text(2*nvar,:),'FontSize',12,'FontWeight','bold')
    
end

axis square
axes=axis;
line([axes(1) axes(2)],[0 0],'LineStyle','--','LineWidth',2,'Color','k')
line([0 0],[axes(3) axes(4)],'LineStyle','--','LineWidth',2,'Color','k')
legend(handle(1:max(class)),clabels,'Location','Best')
legend boxoff

%% Representing pseudo-sample trajectories %%

list=input('For which original variables do you want to represent the corresponding pseudo-sample trajectories? ([1, 2, ...]) ');
pcid=input('For which principal component? ');

figure
hold on

for nvar=1:length(list)
    
    plot(1:np,squeeze(tps_c(:,pcid,list(nvar))),'-')
    text(-0.5,tps_c(1,pcid,list(nvar)),[num2str(list(nvar)),'-'],'Color',cmap_text(2*nvar,:),'FontSize',12,'FontWeight','bold')
    text(np+0.5,tps_c(end,pcid,list(nvar)),[num2str(list(nvar)),'+'],'Color',cmap_text(2*nvar,:),'FontSize',12,'FontWeight','bold')
    
end

axis tight
axes=axis;
axis([-4 np+5 axes(3)-.25 axes(4)+.25])
xlabel('variable range (a.u.)','FontSize',16,'FontWeight','bold')
ylabel(['pseudo-sample score on PC #',num2str(pcid)'],'FontSize',16,'FontWeight','bold')
box off
set(gca,'FontSize',16,'FontWeight','bold','XTick',[1 np],'XTickLabels',{'minimum';'maximum'})
set(gcf,'Color','w')