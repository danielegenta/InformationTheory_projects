% ################
% PART 2 - Exercise 2
% ################

% probability p and vecotr y as input arguments
% Fixed parameters
% alphabet: A = [0 1]
% k = 2

clear all
clc

% GUI
prompt = {'Enter probability p (format: x.x):','Enter vector y (format: x1 x2 x3 x4):'};
dlgtitle = 'Input arguments';
dims = [1 35];
definput = {'0.1','1 1 0 1'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

p = str2double(answer{1});
y = str2num(answer{2});

% may optimize the parameters using I, R matrices (see notebook)
% optimize using bits de2bi(0:card) where card = 2^k
% note: we're using BSC!

G = [1 0 1 1; 0 1 0 1]; % generator matrix (fixed)
% H from G
H = [G(1, 3:4); G(2, 3:4); 1 0; 0 1]; % TO CHECK!

%H = [1 1; 1 0; 1 0; 0 1]; % TO BE GENERATED DINAMICALLY
%y = [1 1 0 1]; % MUST BE AN INPUT - ok
%p = 0.1; % MUST BE AN INPUT - ok 


% #################### PART 1 - BRUTE FORCE APPROACH #######################

% generate codebook 
codebook = [];
for i=0:1
    for j=0:1
        v = [i j]; % generic vector of k bits
        c = v * G; % generic codeword of 2^k bits
        codeword = correctManually(c);
        codebook = [codebook; codeword];
    end
end



% objective (step0): P[c_i = 0 | Y_), P[c_i = 1 | Y_]
    % step1: P[C_ | Y_ ] for each codeword
    % step2: P[c_i | y_i]

% approaching step 1,2
arr_prob_step_2 = [];
for index_cw = 1:length(codebook)

    codeword = codebook(index_cw,:);

    prob_step2 = 1; % P(C_ | Y_ )

    % BSC implementation
    for index_bit = 1:length(codeword)
        if (codeword(index_bit) == y(index_bit))
            prob_step2 = prob_step2 * (1 - p);
        else
            prob_step2 = prob_step2 * (p);
        end
    end  
    arr_prob_step_2 = [arr_prob_step_2 prob_step2];
end

% normalization prob_step2_i/sum(prob_step2_i)
s = sum(arr_prob_step_2);
norm_arr_prob_step_2 = arr_prob_step_2 ./ s; %OK, right results

% approaching step 0
%
% here i compute the final: P[C_i = 0 | Y_] exploiting marginalization
 %P[c_i | Y_]
arr_prob_step_1_0 = [];
arr_prob_step_1_1 = [];
arr_prob_step_1_complete = [];
for j = 1:4
    prob_step_1_0 = 0;
    prob_step_1_1 = 0;
    for i = 1:length(codebook)
        codeword = codebook(i,:);
        prob = norm_arr_prob_step_2(i);
        bit = codeword(j);
        if (bit == 0)
            prob_step_1_0 = prob_step_1_0 + prob;
        end
        if (bit == 1)
            prob_step_1_1 = prob_step_1_1 + prob;
        end
    end
    arr_prob_step_1_0 = [arr_prob_step_1_0 prob_step_1_0];
    arr_prob_step_1_1 = [arr_prob_step_1_1 prob_step_1_1];
end
arr_prob_step_1_complete = [arr_prob_step_1_0; arr_prob_step_1_1];
disp('BRUTE FORCE - probabilities c_i = 0,1 given y (0;1):')
disp(arr_prob_step_1_complete)


% ############################# PART 2 - FACTOR GRAPH APPROACH #################################
disp('##### PART 2 ######');

% step1: given G,find H and translate H into the factor graph (in particular edges)
[rownum,colnum] = size(H);

% manual version
%factor_graph = graph({'F1' 'F1' 'F1' 'F2' 'F2'}, {'C1', 'C2', 'C3', 'C1', 'C4'});
factor_graph = graph();

% ########### GRAPH CONSTRUCTION #########
% add  nodes
% add factors
for i = 1:colnum
    factor_graph = addnode(factor_graph, strcat('F',int2str(i)));
end
% add variables
for i = 1:rownum
    factor_graph = addnode(factor_graph, strcat('C',int2str(i)));
end

%add edges
for row = 1:rownum
    for col = 1:colnum
        if (H(row,col) == 1)
            factor_graph = addedge(factor_graph, {strcat('C',int2str(row))}, {strcat('F',int2str(col))});
        end
    end
end

% H from G
H_bis = [G(1, 3:4); G(2, 3:4); 1 0; 0 1];

% STEP 0
% Compute initial structures
structure_1 = []; % representing c_1 = 0
structure_2 = []; % representing c_1 = 1

% Populating the two structures - may be smarter
for row = 1:rownum
    tuple_str1 = [];
    tuple_str2 = [];
    for col = 1:colnum
        if (H_bis(row,col) == 1)
            if (y(row) == 1)
                tuple_str1 = [tuple_str1 p];
                tuple_str2 = [tuple_str2 1-p];
            else
                tuple_str1 = [tuple_str1 1-p];
                tuple_str2 = [tuple_str2 p];
            end
        else
           tuple_str1 = [tuple_str1 0];
           tuple_str2 = [tuple_str2 0];
        end
    end
    if (row == 1)
        structure_1 = [structure_1 tuple_str1];
        structure_2 = [structure_2 tuple_str2];
    else
        structure_1 = [structure_1; tuple_str1];
        structure_2 = [structure_2; tuple_str2];
    end
end

% STEP 1 - factor graph update
% 1.1 each node sends to the factor p1(0), p1(1) => i exploit structure_1
%                                                   and structure_2

% 1.2 each factor node sends back to the variable 2 messages (see photo)

% 1.3  normalization

iterations = 2;

for ind_iterations = 1:iterations
    structure_2_copy = structure_2;
    structure_1_copy = structure_1;

   for index_factor = 1:col
        for node = 1:row
            if (structure_1(node, index_factor) ~= 0) % node connected to this factor
                %expression
                m = 1;
                for other_node = 1:row
                    if node ~= other_node && structure_2(other_node,index_factor) ~= 0
                        m = m * (1 - 2 * (structure_2(other_node, index_factor)));
                    end
                end

                % espression +1, expression/2
                m = m + 1;
                m = m / 2;

                % update structure
                structure_1_copy(node,index_factor) = m;
                structure_2_copy(node,index_factor) = 1-m;
            end
        end
   end

    % update data structures
    structure_1 = structure_1_copy;
    structure_2 = structure_2_copy;

    % STEP 2 => 2
    % divided in 2 sub-steps: see photo


    % STEP 2.1 => MAP! OUTPUT OF THE PROGRAM (USELESS IN THE FIRST ITERATION)
    % for each variable node in the graph (c1, c2, c3, c4)
    % compute the MARGINAL probability, use all the incoming messages
    % into c_i including what i receive from the channel (GAMMA = 0.9 or 0.1)

    gamma_vector = [p 1-p]; % vector rapresentinginformation from the channel, DO NOT USE P!!!

    a_vec = [];
    b_vec = [];
    % let's start from c_i = 1 => structure 1, GAMMA = 0.1
    for node = 1:row
        a = 1;
        b = 1;
        for factor = 1:col
            if structure_1(node, factor) ~= 0 && structure_2(node, factor) ~= 0 %incoming message into the node_i from factor_j
                a = a * structure_1(node, factor);
                b = b * structure_2(node, factor);
            end
        end
        % marginal probability
        if y(node) == 1
            gamma_a = gamma_vector(1);
        else
            gamma_a = gamma_vector(2);  
        end
        gamma_b = 1-gamma_a;
        a = a * gamma_a;
        b = b * gamma_b;

        % marginal normalized probability
        a_norm = a/(a+b);
        b_norm = b/(a+b);

        a_vec = [a_vec a_norm]; %incoming messages in node c_i, when c_i = 0
        b_vec = [b_vec b_norm];
    end

    if ind_iterations == 2 % END
        break;
    end

    % STEP 2.2
    % each node variable SENDS BACK to the facors the extrinsic information
    for node = 1:row
        for factor = 1:col
            p1 = 1;
            p2 = 1;
            if structure_1(node, factor) ~= 0 % factor and node are linked 
                for factor_star = 1:col
                     if (factor  ~= factor_star) %extrinsic
                         if structure_1(node, factor_star) ~= 0 && structure_2(node, factor_star) ~= 0
                             p1 = p1 * structure_1(node, factor_star);
                             p2 = p2 * structure_2(node, factor_star);
                         end
                     end
                end
                if y(node) == 1
                    gamma_p1 = gamma_vector(1);
                    gamma_p2 = 1-gamma_p1;
                end
                if y(node) == 0
                    gamma_p1 = gamma_vector(2);
                    gamma_p2 = 1-gamma_p1;
                end
                p1 = gamma_p1 * p1;
                p2 = gamma_p2 * p2;

                new_p1 = p1/(p1+p2);
                new_p2 = p2/(p1+p2);

                structure_1_copy(node,factor) = new_p1;
                structure_2_copy(node,factor) = new_p2;
            end       
       end
    end

    structure_1 = structure_1_copy;
    structure_2 = structure_2_copy;
end

output = [a_vec; b_vec]; % row1: c = 0, row2: c=1
disp('DISTRIBUTED COMPUTING - probabilities c_i = 0,1 given y (0;1):')
disp(output)

% since i want binary sum i gotta manually correct (alternatives: mapping)
% (see notebook)
function c = correctManually(c)
    for index=1:length(c)
        if c(index) > 1
            c(index) = 0;
        end
    end
end



