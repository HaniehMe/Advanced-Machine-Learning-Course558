%%%%%%%%%%%%%%%%%%%%%%%
%                                                      %
%            Genetic Algorithm  Code                   %
%                                                      %
%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

%%
tic;
%%
warning('off');

%% Definition of Data Display and Saving parameters

number_x=15; % Maximum x-axis size in the convergence plot of the cost function
S_R=10; % Save results by iteration number
%% paramters setting

CostFunction=@(x,iter) My_Cost(x,iter);       % Cost Function

nvar=3;                    % Number of Variables
xmin= [     0         0         0 ];              % Lower Bound of Variables
xmax=[   500  500   500];               % Upper Bound of Variables

% Modifying the size of xmin and xmax to have a general form
if numel(xmin)==1
  xmin=repmat(xmin,1,nvar);
  xmax=repmat(xmax,1,nvar);
end

%% GA Parameters

npop=25;         % number of population
maxiter=100000;      % max of iteration


pc=0.5;                  % percent of crossover
ncross=2*round(npop*pc/2);   % number of cross over offspring

pm=1-pc;                 % percent of mutation
nmut=round(npop*pm);     % number of mutation offsprig

lb=xmin.*ones(1,nvar); 
ub=xmax.*ones(1,nvar); 

%% initialization

empty.par=[];
empty.fit=[];

elapsed_time=zeros(1,maxiter);

pop=repmat(empty,npop,1);

BESTSolution=nan(maxiter,nvar);
BEST=nan(maxiter,1);
MEAN=zeros(maxiter,1);

%% Display Setting

h_Axes=axes;
h_Plot=plot(h_Axes,1:maxiter,BEST,'LineWidth',2,'YDataSource','BEST');
xlabel(h_Axes,'Iteration');
ylabel(h_Axes,'Cost function Value');
grid(h_Axes,'on');

%%
for i=1:npop
    
   pop(i).par=xmin+rand(1,nvar).*(xmax-xmin);
   pop(i).fit=CostFunction(pop(i).par,i);
   
end

%% main loop

for iter=1:maxiter


% crossover
   
crosspop=repmat(empty,ncross,1);
   
f=sort([pop.fit]);
f=f./sum(f);
f=cumsum(f);

for n=1:2:ncross
   
i1=find(rand<=f,1,'first');
i2=find(rand<=f,1,'first');
        
p1=pop(i1).par;
p2=pop(i2).par;

R=rand(size(p1));

o1=(p1.*R)+(p2.*(1-R));
o2=(p2.*R)+(p1.*(1-R));

crosspop(n).par=o1;                

crosspop(n).fit=CostFunction(o1,n);

crosspop(n+1).par=o2;
            
crosspop(n+1).fit=CostFunction(o2,n);

end
   
% mutation
mutpop=repmat(empty,nmut,1);
  
for n=1:nmut


i=randi([1 npop]); 
    
sol=pop(i).par;

j=randi([1 nvar]);

d=0.1*unifrnd(-1,1)*(xmax(j)-xmin(j));

sol(j)=sol(j)+d;

sol=min(sol,xmax);
sol=max(sol,xmin);

mutpop(n).par=sol;
     
mutpop(n).fit=CostFunction(sol,n);

end
 
% merged
[pop]=[pop;crosspop;mutpop];

% select
[value,index]=sort([pop.fit]);
pop=pop(index);
pop=pop(1:npop);

gpop=pop(1);   % global pop

BESTSolution(iter,:)=gpop.par;
BEST(iter)=gpop.fit;
MEAN(iter)=mean([pop.fit]);

 %%
   if  ( rem(iter,10) == 0) 
      save BESTSolution BESTSolution
      save BEST BEST
      save elapsed_time elapsed_time
   end
    
     %% Display setting
    
   refreshdata(h_Plot);
   drawnow;

    if  iter <= number_x
         axis ( [1, number_x, -inf, inf] );
         set(gca,'XTick',1:number_x)  ;    
     else 
         axis ( [iter-number_x+1, iter, -inf, BEST(iter-number_x+1) ]);
         set(gca,'XTick',iter-number_x+1:iter) ;      
      
    end
    
elapsed_time(1,iter) =toc;
Hour=fix(elapsed_time(1,iter) /3600);
Minute=fix( (elapsed_time(1,iter)-Hour*3600)/ 60 );
Second=fix(elapsed_time(1,iter)-Hour*3600 - Minute*60);
   
title(h_Axes,sprintf('It = %d     BC = %.4e     Time = %d : %d : %d ',iter,BEST(iter),Hour,Minute,Second))  
     
pause(0.01);

end








