clear all
close all
clc
load('iddata-19.mat')

Ts=0.045;
Uid=id.u;
Yid=id.y;
Uval=val.u;
Yval=val.y;
namax=3;
nbmax=3;
nk=1;
m=3;
Nval=length(Uval);
Nid=length(Uid);
y_predictie=zeros(Nval,1);
y_simval=zeros(Nval,1);
y_simid=zeros(Nid,1);
yid_pred=zeros(Nid,1);

na_valori = [];
nb_valori = [];
grad_valori = [];
erori_val = [];
erori_id = [];
mse_heatmap = NaN * ones(namax, nbmax, m);
mse_valori = [];

for grad=1:m
    for na=1:namax
        for nb=1:nbmax
            combinatii=generare_combinatii(na,nb,grad);
            nrcombinatii=length(combinatii);
            d=zeros;
            fi=zeros(Nid,na+nb);
            for k = 1:Nid
                for i=1:na
                    if k-i>0
                        d(k,i)=-Yid(k-i);
                    else
                        d(k,i)=0;
                    end
                end
                for i=1:nb
                    if k-i-nk+1>0
                        d(k,na+i)=Uid(k-i-nk+1);
                    else
                        d(k,na+i)=0;
                    end
                end
                for cont = 1:nrcombinatii
                    fi(k, cont)= 1;
                    for element = 1:na+nb
                         fi(k, cont)= fi(k, cont)*(d(k, element) ^ combinatii(cont, element));
                    end                    
                end
            end

            theta =fi\Yid;
            yid_pred=fi *theta;

            fi_simid = zeros(length(Yid),na+nb);
            vector_id=zeros(Nid, na + nb); 
            for k = 1:Nval
                for i = 1:na
                    if k-i > 0
                     vector_id(k,i)=-y_simid(k-i);
                    end
                end
                for j = 1:nb
                    if k-j-nk+1 > 0
                     vector_id(k,na+j) = Uid(k-j-nk+1);
                    end
                end
                for p = 1 :nrcombinatii
                    fi_simid(k, p) = 1;
                    for el = 1 : na + nb
                        fi_simid(k, p) = fi_simid(k, p) * (vector_id(k, el) ^ combinatii(p, el));
                    end
                end
                y_simid(k) = fi_simid(k,:)*theta;
            end

            fi_simval = zeros(length(Yval),na+nb);
            vector_val=zeros(Nid, na + nb);
            for k = 1:Nval
                for i = 1:na
                    if k-i > 0
                     vector_val(k,i)=-y_simval(k-i);
                    end
                end
                for j = 1:nb
                    if k-j-nk+1 > 0
                     vector_val(k,na+j) = Uval(k-j-nk+1);
                    end
                end
                for p = 1 :nrcombinatii
                    fi_simval(k, p) = 1;
                    for el = 1 : na + nb
                        fi_simval(k, p) = fi_simval(k, p) * (vector_val(k, el) ^ combinatii(p, el));
                    end
                end
                y_simval(k) = fi_simval(k,:)*theta;
            end
            fival = zeros(Nval, nrcombinatii);  
            for k = 1:Nval
                for i=1:na
                    if k-i>0
                        d(k,i)=-Yval(k-i);
                    else
                        d(k,i)=0;
                    end
                end
                for i=1:nb
                    if k-i-nk+1>0
                        d(k,i+na)=Uval(k-i-nk+1);
                    else
                        d(k,i+na)=0;
                    end
                end
                for cont = 1:nrcombinatii
                    fival(k, cont) = 1;
                    for element = 1:na+nb
                        fival(k, cont) = fival(k, cont) * (d(k, element) ^ combinatii(cont, element));
                    end
                end
            end
            y_predictie=fival*theta;

            eroare_valsimulare = mean((Yval(max(na, nb):end) - y_simval(max(na, nb):end)).^2);
            eroare_idsimulare = mean((Yid(max(na, nb):end) - y_simid(max(na, nb):end)).^2);
            eroare_valpredictie = mean((Yval(max(na, nb):end) - y_predictie(max(na, nb):end)).^2);
            eroare_idpredictie = mean((Yid(max(na, nb):end) - yid_pred(max(na, nb):end)).^2);

            mse_valori = [mse_valori; na, nb, grad, eroare_valsimulare, eroare_idsimulare];
            mse_heatmap(na, nb, grad) = eroare_valsimulare;
            na_valori = [na_valori; na];
            nb_valori = [nb_valori; nb];
            grad_valori = [grad_valori; grad];
            erori_val = [erori_val; eroare_valsimulare];
            erori_id = [erori_id; eroare_idsimulare];

            if na==2 && nb==2 && grad==3
                figure
                plot(Yid, 'b'); hold on;
                plot(y_simid, 'r');hold on;
                plot(yid_pred,'g');
                title('identificare')
                legend(['Yid (na=' num2str(na) ', nb=' num2str(nb) ', m=' num2str(grad) ')']);
                figure
                plot(Yval, 'b'); hold on;
                plot(y_simval, 'r');hold on;
                plot(y_predictie,'g');
                title('validare')
                legend(['Yval (na=' num2str(na) ', nb=' num2str(nb) ', m=' num2str(grad) ')']);
            end
       end
   end
end
%%

index = (na_valori == 2) & (nb_valori == 2);
m_fixat = grad_valori(index);
erori_val_fixate = erori_val(index);
erori_id_fixate = erori_id(index);
figure;
plot(m_fixat, erori_val_fixate, '-o', 'LineWidth', 2, 'MarkerSize', 6);
title('Eroare Validare în funcție de m (na = nb = 3)');
xlabel('Gradul m');
ylabel('Eroare Validare');
grid on;
figure;
plot(m_fixat, erori_id_fixate, '-o', 'LineWidth', 2, 'MarkerSize', 6);
title('Eroare Identificare în funcție de m (na = nb = 3)');
xlabel('Gradul m');
ylabel('Eroare Identificare');
grid on;
for grad = 1:m
    figure;
    heatmap(1:namax, 1:nbmax, squeeze(mse_heatmap(:,:,grad)), 'Title', ['MSE pentru gradul m = ', num2str(grad)], ...
        'XLabel', 'na', 'YLabel', 'nb', 'ColorbarVisible', 'on');
end

tabel_mse = array2table(mse_valori, 'VariableNames', {'na', 'nb', 'grad', 'MSE_val', 'MSE_id'});
disp(tabel_mse);

na_unic = unique(na_valori);
nb_unic = unique(nb_valori);

eroare_val_matrice = NaN(length(na_unic), length(nb_unic)); 
eroare_id_matrice = NaN(length(na_unic), length(nb_unic));

for i = 1:length(na_unic)
    for j = 1:length(nb_unic)
        indiceprim = (na_valori == na_unic(i)) & (nb_valori == nb_unic(j)) & (grad_valori == 3);

        if any(indiceprim) 
            eroare_val_matrice(i, j) = mean(erori_val(indiceprim)); 
            eroare_id_matrice(i, j) = mean(erori_id(indiceprim)); 
        end
    end
end

figure
subplot(1, 2, 1); 
[Na, Nb] = meshgrid(nb_unic, na_unic); 
surf(Na, Nb, eroare_val_matrice, 'EdgeColor', 'none'); 
xlabel('Valoare nb');
ylabel('Valoare na');
zlabel('Eroare Validare');
colorbar;
grid on;

subplot(1, 2, 2); 
surf(Na, Nb, eroare_id_matrice, 'EdgeColor', 'none');
xlabel('Valoare nb');
ylabel('Valoare na');
zlabel('Eroare Identificare');
colorbar;
grid on;


function combinatii = generare_combinatii(na, nb, m)
    dimensiuni = na + nb;  
    combinatii = []; 
    function combinatii_recursive(curent, nivel)
        if nivel > dimensiuni 
            if sum(curent) <= m   
                combinatii = [combinatii; curent]; 
            end
            return;
        end
        for valori = 0:m
            combinatii_recursive([curent, valori], nivel + 1);
        end
    end
    combinatii_recursive([], 1); 
end
