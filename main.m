clear all;
close all;

R=10;
samples=1e5;
Theta=1;
Pw=1; %in Watts
Pt_fixed=0.3; %in Watts
Pt_vector= 0:0.1:1; % in Watts
sigma=-5; %in dB
sigma_W=10^(sigma/10);
alpha_vector=0:0.1:5; %in Watts

h = complex(sqrt(0.5)*randn(6, samples),sqrt(0.5)*randn(6, samples));
h = abs(h).^2;

OutageAC=0;
OutageCB=0;

for n = 1: length(Pt_vector)
    Pt=Pt_vector(n);
    
    SIR_AC=((Pt*h(1,:))./(Pw*h(2,:)+sigma_W));
    SIR_CB=((Pt*h(3,:))./(Pw*h(4,:)+sigma_W));
    
    Outage=sum(SIR_AC<Theta | SIR_CB<Theta);
    
    top(n)=(Outage)/samples;
    
    top_theoretical(n)=1-exp(-((Theta*(2*sigma_W))/(Pt)))*(Pt/(Theta*Pw+Pt))^2;
end

for n = 1: length(alpha_vector)
    alpha=alpha_vector(n);
    
    pFA=sum((Pw*h(5,:)+sigma_W)>=alpha);
    pMD=sum((Pt_fixed*h(6,:)+Pw*h(5,:)+sigma_W)<alpha);
    
    dep(n)=(pFA+pMD)/samples;
    
    if(alpha<sigma_W)
        dep_theoretical(n)=1;
    else
        dep_theoretical(n)=1+(Pt_fixed/(Pt_fixed-Pw))*(exp((sigma_W-alpha)/Pw)-exp((sigma_W-alpha)/Pt_fixed));
    end
    
end



%%
% plot the results:
figure(1);
semilogy(Pt_vector, top_theoretical, '-', 'Color', [0.8500 0.3250 0.0980], 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'LineWidth', 1.5, 'MarkerSize',7);
hold on;
semilogy(Pt_vector, top, 's', 'Color', [0.8500 0.3250 0.0980], 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'LineWidth', 1.5, 'MarkerSize',7);
grid on;
ylabel('$\mathrm{Transmission}\,\,\mathrm{Outage}\,\,\mathrm{Probability}$','Interpreter','LaTeX','Fontsize',14);
xlabel('$\mathrm{Covert}\,\mathrm{Transmission}\,\mathrm{Power}$','Interpreter','LaTeX','Fontsize',14);


% plot the results:
figure(2);
semilogy(alpha_vector, dep_theoretical, '-', 'Color', [0.8500 0.3250 0.0980], 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'LineWidth', 1.5, 'MarkerSize',7);
hold on;
semilogy(alpha_vector, dep, 's', 'Color', [0.8500 0.3250 0.0980], 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'LineWidth', 1.5, 'MarkerSize',7);
grid on;
ylabel('$\mathrm{Detection}\,\,\mathrm{Error}\,\,\mathrm{Probability}$','Interpreter','LaTeX','Fontsize',14);
xlabel('$\mathrm{Detection}\,\mathrm{Threshold}$','Interpreter','LaTeX','Fontsize',14);