function feas=OJAG_Wentzell_No_Com(a,K,k,sigma,L,delta,Delta_W_NC)
%% Function which checks the feasibility of the LMIs...
% as derived in Theorem1 i.e with wentzell boundary and leaders communication.

%  a, K,k, sigma are the control gain.
% L is the Lipschitz constant.
% delta is the decay rate.
% Delta_W_C: refers to the spacing between two consecutive leaders. 

% Decision variables and notations 
sdpvar mu lambda_1  lambda_2 lambda_3 lambda_4 p1 p2 p3

% Matrix P > 0. 
P=blkvar;
P(1,1) =p1;
P(1,2) = p2; 
P(2,2)= p3;
P=sdpvar(P); 

% Matrix N < 0. 
N=blkvar;
N(1,1) =-(2*(K-delta)-lambda_2*L^2);
N(1,4) = 1; 
N(1,5)= K;
N(2,2) =- (2*mu*(K-delta)+2*a-lambda_1*((Delta_W_NC)^2/pi^2)) ;
N(3,3) = - 2*mu*a; 
N(3,4)= -mu; 
N(3,5) =  -mu*K; 
N(4,4) = -lambda_2;
N(5,5) = -lambda_1;
N=sdpvar(N); 

% Matrix M < 0. 
M=blkvar;
M(1,1) = -2*(k +K- delta)*p1  + lambda_3*L^2;
M(1,2) = -2*(k +K- delta)*p2;
M(1,3) =  sigma*p1 + mu*k -a ; 
M(1,4) = -sigma*p2;
M(1,5) = p1;
M(1,6)= p2;
M(2,2) = -2*(k + K - delta)*p3 + lambda_4*L^2;
M(2,3) = sigma*p2;
M(2,4) = -sigma*p3-mu*k+a  ;
M(2,5)= p2;
M(2,6) = p3; 
M(3,3) = -2*sigma*mu;
M(3,5) = -mu;
M(4,4) = -2*sigma*mu;
M(4,6)= mu;
M(5,5) = -lambda_3;
M(6,6) = -lambda_4;
M=sdpvar(M); 


% Solution of LMIs
LMIs=[mu >=0,lambda_1>=0,lambda_2>=0, lambda_3>=0, lambda_4>=0,P>=0, N<=0, M<=0 ];
options=sdpsettings('solver','sedumi','verbose',0);
sol=optimize(LMIs,[],options); 

feas=0; 
if sol.problem==0
    [primal,~]=check(LMIs); 
    feas=min(primal)>=0; 
else
    yalmiperror(sol.problem) 
end