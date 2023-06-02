% clc;close all; 
load('ALC_DL_DO_150.mat','C')    % Load the codebook
M=size(C,2);  % size of the codebook for one user
m=log2(M); % no of bits in a symbol/block
K=4;   % No of orthogonal PREs/ OFDM subcarriers
J=6;   % No of users/layers
Nr=4; % No of BS antenaas
F=get_indicator_matrix_factor_graph(C, J, K);
d_v=length(find(F(:,1)));
d_f=length(find(F(1,:)));
%% %%%%%%% Power Allocation %%%%%%%%%
power_users=ones(1,J);  
sqrt_power_matrix=ones(K,J);
%% Monte Carlo %%%%%%%%%%%%%%%%%%%%%%%%
%Eb_N0_dB=30;
Eb_N0_dB=20;
%Eb_N0_dB=5;
%Eb_N0_dB=20;
Eb_N0=10.^(Eb_N0_dB/10);
Es=sum(sum((abs(C)).^2))/length(C);
Eb=Es/m;   
N0=Eb./Eb_N0;
SNR=Es./N0;
SNR_db=10*log10(SNR);
sigma=sqrt(N0/2);%standard deviation of real gaussian
%max_block_errors=[  3000 2500 2300 2100  1700 1600   ];
max_block_errors=  10;
%max_block_errors= [10 10 10 10 10 10 10 ]   ;

%al=[0.3, 0.2, 0.4,0.3,0.4,0.7, 0.4,0.6,0.7,0.6];
%al=[45,70,75,80,90,100];
T=30;
al=300;
%rho=5;
%CEE=0:0.1:0.5; 
SER=zeros(1,length(al));
%% 
 for e=1:length(al)
    %rho=N0/Es(e);
    rho=N0/Es;
    block_errors=0;
    symbol_errors=0;
    block=0;
    while block_errors<max_block_errors(e)
        %%   SCMA Encoding %%%%%%
        bits=randi([0 1],J,m);% blocks of bits for all users
        %mapping
        symbol_indices=bi2de(bits,'left-msb')+1; % symbols for all users
        SCMA_codewords=zeros(K,J);  % collection of the codewords for all users
        al_R=zeros(1,J);
        be_I=zeros(1,J);
        for j=1:J         % for each user
            present_codebook=C((j-1)*K+1:j*K,:);   % codebook for the jth user
            SCMA_codewords(:,j)=present_codebook(:,symbol_indices(j));
            al_R(1,j)=max(max(abs(real(present_codebook))));
            be_I(1,j)=max(max(abs(imag(present_codebook))));
        end
        
        %% Transmission through Rayleigh fading channel %%
        AWGN=sigma*(randn(Nr*K,1)+1j*randn(Nr*K,1));  
        % complex Gaussian noise
        %AWGN=0;
        H=zeros(Nr*K,J*d_v);
        %H_cap=zeros(Nr*K,J*d_v);
       SCMA_CW_MU=zeros(J*d_v,1);
       %e1=0;
       for j=1:J
            h=zeros(Nr*K,d_v);
            %h_cap=zeros(Nr*K,d_v);
            a=find(~F(:,j));%zero element indices
            for nr=1:Nr
                h_j=diag(1/sqrt(2)*(randn(K,1)+1j*randn(K,1)));                
               
                %h_j_cap=h_j+CEE(e)*diag(1/sqrt(2)*(randn(K,1)+1j*randn(K,1)));
                
                h_j(:,[a(1) a(2)])=[];% 150% overloading
                %h_j_cap(:,[a(1) a(2)])=[];
                
                h((K*nr-(K-1)):K*nr,:)=h_j;% 150% overloading
                %h_cap((K*nr-(K-1)):K*nr,:)=h_j_cap;
               
            end
            H(:,d_v*j-1:d_v*j)=h;
            %H_cap(:,d_v*j-1:d_v*j)=h_cap;
            SCMA_CW=SCMA_codewords(:,j);%j-th user codeword
            SCMA_CW([a(1) a(2) ])=[];% 150% overloading
           
            SCMA_CW_MU(d_v*j-1:d_v*j,1)=SCMA_CW;           
       end       
        y=H*SCMA_CW_MU+AWGN;
       % y=H*SCMA_CW_MU;
        %% received SCMA codeword UP_link               
        [SCMA_MU_CW_det]=DPS_ADMM_SCMA_150(H,y,J,rho,al(e),T,d_v,al_R,be_I);%150% overloading
        %SCMA_MU_CW_det=MMSE_SCMA_150(H,y,J,rho,d_v);
       
        det_ind=zeros(J,1);
        for j=1:J
            SCMA_CW_det=SCMA_MU_CW_det(2*j-1:2*j);
             present_codebook=C((j-1)*K+1:j*K,:); %j-th user codebook
             a=find(~F(:,j));
             present_codebook([a(1) a(2)],:)=[];% 150% overloading
            
             ED1=zeros(1,M);
             for m1=1:M
                 ED=norm(SCMA_CW_det- present_codebook(:,m1));
                 ED1(1,m1)=ED;
             end
           det_ind(j,1)=find(ED1==min(ED1));           
                
        end
        
        error_locations=find(symbol_indices~=det_ind);
        if ~isempty(error_locations)
            block_errors=block_errors+1;           
            symbol_errors=symbol_errors+length(error_locations);
            %fprintf('\n  %d error collected',block_errors); 
        end     
        block=block+1;
    end
    SER(e)=symbol_errors/block/J;% Each block contains J-symbols
    %BER(e)=SER(e)/m;
    fprintf('\n Simulation done for %d alpha',al); 
end
semilogy(al,SER,'b-*','LineWidth',2) ;
xlabel('al');
ylabel('SER');
legend('al Vs SER');





