function [recovered_data] = MF_receiver(noisy_message, NumPts, h_opt_s1)
    
    recovered = filter (h_opt_s1, 1,noisy_message)/NumPts ;
    
    sample_rec = recovered (NumPts:NumPts:end); 
    decode_emg = sign(sample_rec); 
    decode_emg(decode_emg == -1) = 0;
    t_1 = rem(length(decode_emg),7);
    decode_emg = decode_emg (1:end-t_1);
    
    reshape_emg = reshape(decode_emg,[7, length(decode_emg)/7]) ;
    emg_msg = int2str(reshape_emg)';
    emg_msg = emg_msg(1:3:end,:);
    recovered_data = char(bin2dec(emg_msg))';
                                    
end