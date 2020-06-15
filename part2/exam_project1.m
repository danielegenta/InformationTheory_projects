% ################
% PART 2 - Exercise 1
% ################

% STEP1: Generate random sequence of l bits
% output: random l-bit sequence

clear all
clc

ber_arr = [];
ber = 1;
SNR = 0;
input_number_wrong_bits = 1000;
% SNR = O, ..., 6
while SNR < 5
    
    number_wrong_bits = 0;
    l = 0;
    l_considered = 1;
    
    % information bits number depends on the number of wrong bits output of
    % Viterbi
    while number_wrong_bits < input_number_wrong_bits

        l = 100; % parameter, dynamically change until i have at least input_number_wrong_bits wrong bits
        rand_bin = round(rand(1,l));

        % Intermediate STEP: zero-pad l bits sequence with 10 dummy values at the
        % end
        dummy = zeros(1,10);
        final_bit_seq = [rand_bin dummy];

        % STEP2: Convolutional TRELLIS encoder 
        % input: (l+10) bits
        % output: (l+10) * n bits where n is the number of states of the convolutional
        %           TRELLIS encoder
        trellis_input_2states = [[0,0], [0,1]; [1,1], [1,0]]; % trellis diagram corresponding matrix: row_index = state (1 or 2), column index = input, value = array of bits (output)

        output = zeros(1, 2*length(final_bit_seq));
        o_ind = 1;
        state = 0; % Supposition: initially i'm in the zero state (QUESTION)

        for index = 1:length(final_bit_seq)

            % initial state
            input_i = final_bit_seq(index);

            if (index == 1)
                state = 1;
            end

            if input_i == 0
                start = 1;
            else
                start = 3;
            end

            output(o_ind) = trellis_input_2states(state, start);
            output(o_ind + 1) = trellis_input_2states(state, start+1);

             % next state
            if input_i == 0
                state = 1;
            else
                state = 2;
            end

            o_ind = o_ind + 2;

        end

        %STEP3: 2-PSK modulator
        output_psk = [];
        %alpha = -1,+1
        for index = 1:length(output)
          input_i = output(index);
          if input_i == 0
              alpha = -1;
          else
              alpha = 1;
          end
          output_psk = [output_psk alpha];
        end

        %STEP4: AWGN (first implementation)
        % this block adds to each symbol (-alpha, alpha) a sample
        % of a Gaussian rv with zero mean and variance sigma^2 = N_0 /2

        % parameters: N_0 = 1/SNR (since Eb = 1), sigma = sqrt(N_0/2)
        % i expect that with larger SNR, smaller BER

        %      i should have an outer cycle in order to compute the L rand bits
        %      multiple times with different SNR values (from zero to six) 
        %      (STOP CONDITION: at least 100 wrong info bits after viterbi) 
        %      (!!!) from theory: i expect a better BER using a larger SNR

        % SNR is measured in dB in order to be used to compute N_0 it must be
        % converted into a number using the formula: SNR_norm =

        %SNR = 6; from cycle
        SNR_norm = 10^(SNR/10);    %check formula
        N_0 = 1/SNR_norm; % variance
        stdev = sqrt(N_0/2);
        noise_vector = stdev .* randn(1,length(output_psk));
        output_awgn = output_psk + noise_vector;
        disp(noise_vector);

        %STEP5: 2-PSK decision rule
        output_psk_decision = [];
        for index = 1:length(output_awgn)
            i = output_awgn(index);
            dist = abs(1-i);
            o = -1;
            if (dist < 0.5)
                o = 1;
            else
                o = 0;
            end
             output_psk_decision = [output_psk_decision o];
        end

        % STEP6: Viterbi Algorithm
        % input: output_psk_decision, input
        % output: L information bits
        % compare L input bits and L information bits

        % COST per edge - Hamming distance 

        % suggestion: save information bits during Viterbi
        % not dynamic: i imagine to have just 2 states

        state_metric_trajectories = [0 0]; %[state_metric_traj_1 state_metric_traj_2]

        info_bits_traj_1 = [];
        info_bits_traj_2 = [];

        index_bits = 1;
        initial_state = 1;

        state_trajectories = [-1 -1]; % [state_traj_1 state_traj_2]

        for index = 1:length(final_bit_seq)

            % fetch the right 2 bits from the viterbi input sequence
            bits_considered = [output_psk_decision(index_bits) output_psk_decision(index_bits+1)];
            index_bits = index_bits + 2; 

            % initial state => 2 trajectories created
            if (index == 1)
                % traj.1
                edge_considered_t1 = [trellis_input_2states(initial_state, 1) trellis_input_2states(initial_state, 2)];
                cost1 = hammingDistance(bits_considered, edge_considered_t1);
                state_metric_trajectories(1) = state_metric_trajectories(1) + cost1;
                state_trajectories(1) = 1;
                info_bits_traj_1 = [info_bits_traj_1 (state_trajectories(1)-1)];

                % traj.2
                edge_considered_t2 = [trellis_input_2states(initial_state, 3) trellis_input_2states(initial_state, 4)];
                cost2 = hammingDistance(bits_considered, edge_considered_t2);
                state_metric_trajectories(2) = state_metric_trajectories(2) + cost2;
                state_trajectories(2) = 2;
                info_bits_traj_2 = [info_bits_traj_2 (state_trajectories(2)-1)];
            end

            if (index > 1)

                % for each state of each trajectory i have 2 possible paths, i
                % choose the one with the minimum cost, if the cost is the same i
                % choose randomly

                for index_state = 1:length(state_trajectories) % iterating over the 2 trajectories
                    state_traj_i = state_trajectories(index_state);

                    % dinamically work either if the state is 0,1 (row 1 or
                    % 2 of trellis)
                    % edge 1 cost computation 
                    first_edge_considered = [trellis_input_2states(state_traj_i, 1) trellis_input_2states(state_traj_i, 2)];
                    first_edge_cost = hammingDistance(bits_considered, first_edge_considered);

                    % edge 2 cost computation 
                    second_edge_considered = [trellis_input_2states(state_traj_i, 3) trellis_input_2states(state_traj_i, 4)];
                    second_edge_cost = hammingDistance(bits_considered, second_edge_considered);

                    % choose edge
                    if (first_edge_cost < second_edge_cost)
                        % update state 
                        state_trajectories(index_state) = 1;
                        % update state metric 
                        state_metric_trajectories(index_state) = state_metric_trajectories(index_state)+first_edge_cost;
                    end
                    if (second_edge_cost < first_edge_cost)
                        state_trajectories(index_state) = 2;
                        state_metric_trajectories(index_state) = state_metric_trajectories(index_state)+second_edge_cost;
                    end
                    % (in case of parity i choose one randomly)
                    if first_edge_cost == second_edge_cost
                        r = rand;
                        if (rand < 0.5)
                            state_trajectories(index_state) = 1;
                            state_metric_trajectories(index_state) = state_metric_trajectories(index_state)+first_edge_cost;
                        else
                             state_trajectories(index_state) = 2;
                            state_metric_trajectories(index_state) = state_metric_trajectories(index_state)+second_edge_cost;
                        end
                    end

                     % update info bits
                    if (index_state == 1) %traj1
                        info_bits_traj_1 = [info_bits_traj_1 (state_trajectories(index_state)-1)];
                    else %traj2
                        info_bits_traj_2 = [info_bits_traj_2 (state_trajectories(index_state)-1)];
                    end
                end 
            end
        end

        % select final winning trajectory
        if (state_metric_trajectories(1) < state_metric_trajectories(2))
            winning = info_bits_traj_1;
        else
            winning = info_bits_traj_2;
        end


        %STEP7: Compute number of wrong bits, BER (#wrong/#rx),
        % SNR, ALPHA fixed.
        % note: do not consider dummy bits (last 10 bits)!

        % i do not consider dummy bits
        number_wrong_bits_i = hammingDistance(winning(1:l), final_bit_seq(1:l));
        number_wrong_bits = number_wrong_bits + number_wrong_bits_i;
        ber = number_wrong_bits/(l_considered*l); % do not consider last 10 dummy bits
        l_considered = l_considered+1;
    end
    SNR = SNR + 1;
    %ber = number_wrong_bits/length(winning(1:l));
    ber_arr = [ber_arr ber];
end
x = linspace(0,5,5);
y = ber_arr;
plot(x,y);
ylabel('BER');
xlabel('SNR');
%disp(ber_arr);

% ########################## UTILITY FUNCTIONS ###########################
% Hamming Distance function
% input:
% output:
% theory bg: compute number of different bits between 2 bit sequences
function distance = hammingDistance(v1, v2)
    % can't use without a specific extension
    %matrix = [v1; v2];
    %distance = size(matrix,2)*pdist(matrix,'hamming');
    distance = 0;
    for index = 1:length(v1)
        if (v1(index) ~= v2(index)) % different bit at position index
            distance = distance + 1;
        end
    end
end
