
% Set diagonal block as the in-network edges and off diagonal as
% out-network
cor_f=csvread('cor_fvicc.csv');

edge_ind_matrix=zeros(184,184);
for i=1:8
    
    edge_ind_matrix(CindxVICC==CIDVICC(i),CindxVICC==CIDVICC(i))=1;
end
for i=1:184
    edge_ind_matrix(i,i)=0;
end
in_edge_ind=find(squareform(edge_ind_matrix)==1);
out_edge_ind=find(squareform(edge_ind_matrix)==0);
CorZVICC_n0_in=cor_f(in_edge_ind);
CorZVICC_n0_out=cor_f(out_edge_ind);


% Plot the distribution of the in and out transformed correlation
histogram(CorZVICC_n0_in,40);
hold on
histogram(CorZVICC_n0_out,40)
legend("Inside-network edges","Outside-network edges")
xlabel("Fisher's Z transformed correlation coefficient")

csvwrite('all.csv',cor_f);
csvwrite('in.csv',CorZVICC_n0_in);
csvwrite('out.csv',CorZVICC_n0_out);

% locfdr_in_n_out.R
edge_length=length(cor_f);
locfdr_all_in_out=csvread('fdr_p0.csv');
p0_all=locfdr_all_in_out(edge_length+1);
p0_in=locfdr_all_in_out(edge_length+2);
p0_out=locfdr_all_in_out(edge_length+3);
locfdr_list=locfdr_all_in_out(1:edge_length);
thres_in=1/(4*p0_in/(1-p0_in)*(1-p0_all)/p0_all+1);
thres_out=1/(4*p0_out/(1-p0_out)*(1-p0_all)/p0_all+1);

locfdr_thres=locfdr_list;
locfdr_thres(in_edge_ind)=locfdr_list(in_edge_ind)<=thres_in;
locfdr_thres(out_edge_ind)=locfdr_list(out_edge_ind)<=thres_out;
VICC_thres=squareform(locfdr_thres);
figure;imagesc(VICC_thres(ClistVICC,ClistVICC));


% Compare with glasso
Cor_node=squareform(CorZ1924_n0');
node_num = size(Cor_node,1);
for i = 1:node_num
    Cor_node(i,i)=1;
end