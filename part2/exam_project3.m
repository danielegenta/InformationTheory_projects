% ################
% PART 2 - Exercise 3
% ################

clear all
clc

% OUTPUT FORMAT
% Set of rules
% each rule has the format [x1_val x2_val x3_val x4_val]

% training_set = [{0,0,10,0};
%                 {0,0,70,0};
%                 {0,1,20,0};
%                 {0,1,80,1};
%                 {1,0,40,0};
%                 {1,0,60,1};
%                 {1,1,50,0};
%                 {1,1,60,1};
%                 ];

% x1, x3 are numerical
training_set = [{30,0,10,0};
                {30,0,70,0};
                {30,1,20,0};
                {30,1,80,1};
                {60,0,40,0};
                {60,0,60,1};
                {60,1,50,0};
                {60,1,60,1};
                ];
            
numberOfFeatures = numel(training_set(1,:));

% dinamically create labels for the training set (useful in tree creation)
% labels = [x1 x2 x3 ... xn class]
labels = cell(1, numberOfFeatures);
for i=1:numberOfFeatures-1
    labels{i} = (strcat('x',num2str(i)));
end
labels{numberOfFeatures} = 'class';
            
% class measures      
class_index = numel(training_set(1,:));
class_vector = training_set(:, class_index);
class_probs = extractProbabilities(class_vector);
class_entropy= H(class_probs);

% i save the used tresholds in order to normalize the input vector by the
% user
tresholds = []; 
    
% -- ID3 IMPROVEMENT --
% step-0: detect numerical features
numberOfFeatures = numel(training_set(1,:)); 
for i = 1:numberOfFeatures-1 % i do not consider the class
    q = cellfun(@(x) isnumeric(x) && numel(x)==1, training_set);   % true for elements of C that are numerical scalars
    IsAllNum = all(q,1);
    
    tmp_training_set = training_set;
    
    % completely numeric vector
    if IsAllNum == 1
        array = cell2mat(training_set(:, i));
        if numel(unique(array)) > 2 
            % not interested in binary columns (?)
            % work on i numeric column (column 3 in the example)
            % fixing the treshold t, the t is 
            % each value found of the numeric column
            cell_array_binarized_column = [];
            igrs = [];
            
            max_igr = 0;
            max_training_set = {};
            max_t = 0;
            for j = 1:numel(array)
                t = array(j);
                binarized_column = binarization(array, t);
                
                % substitute binarized column in the training set
                % compute for each treshold the IGR
                tmp_training_set(:,i) = binarized_column;
                igr = computeSingleIgr(tmp_training_set, class_index, class_entropy, i);
                if (igr > max_igr)
                    max_igr = igr;
                    max_training_set = tmp_training_set;
                    max_t = t;
                end
                %igrs = [igrs igr]; %useful just to debug
                training_set = max_training_set;
                %
                
            end
            tresholds = [tresholds max_t];
            % select which configuration gives me the best IGR
             igrs;
        % useful for the first column in the training set
        else
            % check if the numbers are 0,1 or else binarize the column
            binary_vec = [0,1];
            C=intersect(array,binary_vec);
            if numel(C)<2
                t = min(array);
                binarized_column = binarization(array, t);
                tmp_training_set(:,i) = binarized_column;
                training_set = tmp_training_set;
                tresholds = [tresholds t];
            end
        end
        
    
    end
end
% -- END ID3 IMPROVEMENT --

% from now on i can apply the same computation as ex2
% ------ TRAINING -------
% working on class        
class_index = numel(training_set(1,:));
class_vector = training_set(:, class_index);
class_probs = extractProbabilities(class_vector);
class_entropy= H(class_probs);

rule = {};
setOfRules = {};

% the output of this function is a set of logical rules that map a decision tree classifier,
% each of them has this format:
% [x1 x2 x3 class] = [val1 val2 val3 class_value]
% ex: [null no null yes]
setOfRules = createDecisionTree(training_set, class_index, class_entropy, labels, numberOfFeatures,rule, setOfRules);
% set of logical rules in the format: [x1_val x2_val x3_val class_val]
disp(setOfRules);

% TESTING
% GUI
prompt = {'Enter vector v to classify (format: x1 x2 x3):'};
dlgtitle = 'Input arguments';
dims = [1 35];
definput = {'30 1 20'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
v1 = str2num(answer{1});
normalized_v1 = normalizeInputVector(v1, tresholds);
if numel(normalized_v1) == (size(training_set,2)-1)
    c1 = assignClass(normalized_v1, setOfRules);
    fprintf('predicted class = %d\n', c1{1});
end

function v_new = normalizeInputVector(v1, tresholds)
    v_new = [];
    cnt = 1;
    for i = 1:numel(v1)
        if v1(i) > 1 % i do not binarize 0,1 ..
            if v1(i) <= tresholds(cnt)
                v_new = [v_new 0];
            else
                v_new = [v_new 1];
            end
            cnt = cnt+1;
        else
            v_new = [v_new v1(i)];
        end
    end
end

% function used to binarize a given feature based on a certain treshold
function binarized = binarization(column, t)
    binarized = {};
    for el = 1:numel(column)
        if column(el) <= t
            new = 0;
        else
            new = 1;
        end
        binarized{el} = new;
    end
end

% extract probabilities from vector of categorical attributes
function p = extractProbabilities(v)
    array = cell2mat(v);
    x = unique(array);
    p = zeros(1,length(x));
    for i = 1:length(x)
        p(i) = sum(array==x(i))/length(array);
    end 
end

% calculate entropy
function h = H(p)
    p(p==0) = [];
    h       = -sum(p.*log2(p));
end

function igr_c_x_i = computeSingleIgr(training_set, class_index, class_entropy, feature_i)
    x_i = training_set(:, feature_i);
    x_i_prob = extractProbabilities(x_i);
    x_i_entropy = H(x_i_prob);
    % conditional entropies for each possible attribute of x1
    G = cell2mat(x_i);
    GN = unique(G);
    conditional_entropies_x_i = [];
    for el = 1:numel(GN)
        conditional_training_set = [];
        for j = 1:numel(training_set(:,feature_i))
            if (training_set{j, feature_i} == GN(el))
                conditional_training_set = [conditional_training_set; training_set(j, :)];
            end
        end
        conditional_probs = extractProbabilities(conditional_training_set(:, class_index));
        conditional_entropy = H(conditional_probs);
        conditional_entropies_x_i = [conditional_entropies_x_i conditional_entropy];
    end
    % conditional entropy H(C|X1) as weighted sum of the conditional entropies
    x_i_cardinality = numel(GN);
    w_sum = 0;
    for index = 1:numel(GN)
        % extract weight
        array = cell2mat(training_set(:, feature_i));
        w_i = sum(array==GN(index))/length(array);
        w_entropy = conditional_entropies_x_i(index)*w_i;
        w_sum = w_sum + w_entropy;
    end
    x_i_c_entropy = w_sum;
    % mutual information
    mi_c_x_i = class_entropy - x_i_c_entropy;
    % information gain ratio
    igr_c_x_i = mi_c_x_i/x_i_entropy;
end



% ---- TRAIN FUNCTIONS -------
% train the decision tree
function setOfRules = createDecisionTree(training_set, class_index, class_entropy, labels, numberOfFeatures, rule, setOfRules)

    igrs = Extract_igr(training_set, class_index, class_entropy);
    [argvalue, argmax] = max(igrs);
    % find possible values of the argmax feature
    argmax_column = cell2mat(training_set(:, argmax));
    [GN, ~, G] = unique(argmax_column);
    % iterate over them in order to get new datasets
    % split the dataset according to the selected feature
    feature_examined = str2double(extractBetween(labels{argmax},2,2));
    
    alreadyDeletedLabel = false;
    for i = 1:numel(GN)
        
        % processing GN(i)
        value_examined = GN(i);
        rule{feature_examined} = value_examined;
     
        new_training_set = [];
        for j = 1:numel(training_set(:,1))
            if (training_set{j, argmax} == GN(i))
                new_training_set = [new_training_set; training_set(j, :)];
            end
        end
        % my logical rule: [x2=yes  x1=very_high class=yes]
        
        % i delete the already used feature from the new training sets
        new_training_set(:,argmax) = [];
        
        
        if (~alreadyDeletedLabel)
            labels{1, argmax} = [];
            labels = labels(~cellfun('isempty',labels));
            alreadyDeletedLabel = true;
        end
        class_index = numel(new_training_set(1,:)); %upd

        % STOPPING CRITERION - 1
        % detect leaf nodes
        class_column = cell2mat(new_training_set(:, class_index));
        if all(class_column == class_column(1))
            % thats a leaf node
            %      map into logical associations
            %      deploy the rule
            initial_class_index = numberOfFeatures;
            rule{initial_class_index} = class_column(1); % class
            if numel(setOfRules) == 0
                setOfRules = rule;
            else
                setOfRules(end+1, :) = rule;
            end
            rule{feature_examined} = [];

        else
            % recursive step
            %TODO check if i have just one remaining attribute whose value
            % is not to be splitted
            go_on = true;
            features_cardinality = numel(new_training_set(1,:))-1;
            if (features_cardinality == 1)
                feature_column = cell2mat(new_training_set(:, 1));
                if all(feature_column == feature_column(1))
                    % leaf node
                    go_on = false;
                end
            % STOPPING CRITERION - 2 - TO TEST
            elseif (features_cardinality == 0)
                % detect the most common class among the remaining vectors
                % leaf node labeled by the most common class among the
                % subset vectors
                most_common_class = MostCommonClass(class_column);
                
                initial_class_index = numberOfFeatures + 1;
                rule{initial_class_index} = most_common_class; % class
                if numel(setOfRules) == 0
                    setOfRules = rule;
                else
                    setOfRules(end+1, :) = rule;
                end
                rule{feature_examined} = [];
                
                go_on = false;
            end
            if go_on
                setOfRules = createDecisionTree(new_training_set, class_index, class_entropy, labels, numberOfFeatures,rule, setOfRules);
            end
        end
    end
end

function igr_vector = Extract_igr(training_set, class_index, class_entropy)
    % working on feature i
    % from all features select the one with the larger IGR
    igr_vector = [];
    for feature_i = 1:numel(training_set(1,:))-1
        x_i = training_set(:, feature_i);
        x_i_prob = extractProbabilities(x_i);
        x_i_entropy = H(x_i_prob);
        % conditional entropies for each possible attribute of x1
        G = cell2mat(x_i);
        GN = unique(G);
        conditional_entropies_x_i = [];
        for el = 1:numel(GN)
            conditional_training_set = [];
            for j = 1:numel(training_set(:,feature_i))
                if (training_set{j, feature_i} == GN(el))
                    conditional_training_set = [conditional_training_set; training_set(j, :)];
                end
            end
            conditional_probs = extractProbabilities(conditional_training_set(:, class_index));
            conditional_entropy = H(conditional_probs);
            conditional_entropies_x_i = [conditional_entropies_x_i conditional_entropy];
        end
        % conditional entropy H(C|X1) as weighted sum of the conditional entropies
        x_i_cardinality = numel(GN);
        w_sum = 0;
        for index = 1:numel(GN)
            % extract weight
            array = cell2mat(training_set(:, feature_i));
            w_i = sum(array==GN(index))/length(array);
            w_entropy = conditional_entropies_x_i(index)*w_i;
            w_sum = w_sum + w_entropy;
        end
        x_i_c_entropy = w_sum;
        % mutual information
        mi_c_x_i = class_entropy - x_i_c_entropy;
        % information gain ratio
        igr_c_x_i = mi_c_x_i/x_i_entropy;

        igr_vector = [igr_vector igr_c_x_i];
    end
end

% ------------ TEST FUNCTION ---------------
function c = assignClass(v, setOfRules)
    for rule_idx = 1:numel(setOfRules(:,1))
        rule = setOfRules(rule_idx, :);
        matching = isMatching(v, rule);
        class_idx = numel(rule);
        if (matching)
            c = rule(class_idx);
            break;
        end
    end
end

function bol = isMatching(v, rule)
    bol = true;
    for feature_idx = 1:numel(rule)-1
        if  ~isempty(rule{feature_idx})
            if rule{feature_idx} ~= v(feature_idx)
                bol = false;
                break;
            end
        end
    end
end