clear
clc

%Gutenberg-Richter relationship with 𝑏 = 1
%code to generate the hazard curve

Mmin = 6;
Mmax = 8; 

m = 6:0.1:8;

b =1;

Fm = zeros(size(m,2),1);
for i = 1:size(m,2)
    Fm(i) = (1- 10^(-b*(m(i) - Mmin)))/(1- 10^(-b*(Mmax - Mmin)));
end

Pm = zeros(size(m,2),1);
for i = 1:size(m,2)
    if i == length(m)
        break
    else
        Pm(i) = Fm(i+1) - Fm(i);
    end
end
figure; 
bar(m,Pm)
saveas(gcf,'barchart.png')

fm = zeros(size(m,2),1);
for i = 1:size(m,2)
    fm(i) = (b*log(10)*10^(-b*(m(i) - Mmin)))/(1- 10^(-b*(Mmax - Mmin)));
end
figure;
plot(m, fm)

%CHARACTERISING DISTRIBUTION OF DISTANCES TO FUTURE EARTHQUAKE
%cdf for the distance R
R = 25;
fault_length = 120;
hypoR = sqrt(R^2 + (fault_length/2)^2);
r = 0.5:0.5:fault_length/2 + 10;

for i = 1:size(r,2)
    if r(i) < R
        F(i) = 0;
    elseif r(i) >= R && r(i) < hypoR
        F(i) = (2*sqrt(r(i)^2 - R^2))/fault_length;
    else
        F(i) = 1;
    end
end
figure;
plot(r, F)
saveas(gcf, 'cdf.png')
%pdf for the distance R
for i = 1:size(r,2)
    if r(i) < R
        f(i) = 0;
    elseif r(i) >= R && r(i) < hypoR
        f(i) = r(i)/(fault_length/2 * sqrt(r(i)^2 - R^2));
    else
        f(i) = 0;
    end
end
figure;
plot(r, f)
saveas(gcf, 'pdf.png')

% Predicting Ground Motion Intensity
PGA = 0.01:0.01:6;

e1 = 0.4383; 
e2 = 0.106;
c1 = -0.5543;
c2 = 0.0195; 
c3 = -0.0075; 
gamma = 0.01; 

std_SA = 0.6; 
%mean of the PGA (Sa)
for i = 1:size(m,2)
    for j = 1:size(r,2)
            mean_InSA(i,j) = e1 + e2 *(m(i) - 6.75) + (c1 + c2 *(m(i) - 4.5))*log(r(j)) + c3 *(r(j)-1);
    end
end

%coverting the mean_InSA into a row/column matrix
convert_meanInSA = zeros(size(mean_InSA,1)*size(mean_InSA,2), 1);
u = 1;
for i = 1:size(mean_InSA,1)
    convert_meanInSA(u:u+size(mean_InSA,2)-1, 1) = mean_InSA(i,:);
    u = u+size(mean_InSA,2);

end

%cell_PGA = cell(size(PGA,2),1);
%P(PGA > Xg) for each PGA from 0.01 to 6
for i = 1:size(PGA,2)
    lnPGA = log(PGA(i));
        for k = 1:length(convert_meanInSA)
            cell_PGA(i, 1:k) = 1 - normcdf((lnPGA - convert_meanInSA(k))/std_SA);
        end 
end  

%saving the result of each PGA into a matrix and kept in a cell

for i = 1:size(PGA,2)
    transpose_PGA = cell_PGA';
    cell_PGA2{i} = reshape(transpose_PGA(:,i), [size(mean_InSA,1), size(mean_InSA,2)]);
end


for j = 1:size(cell_PGA,1)
    ProductPGA = Pm .* cell_PGA2{j};
    gammaPGA(j) = gamma*sum(ProductPGA)*F';

end

figure;
semilogy(PGA, gammaPGA)
ylabel('Annual Rate of Exceedance ')
xlabel('PGA (g)')
ylim([gammaPGA(length(gammaPGA)),1e0])

spectraSDBE = interp1(gammaPGA, PGA, 0.1);
spectraSAMCE = interp1(gammaPGA, PGA, 0.02);

disp(spectraSAMCE)
disp(spectraSAMCE)