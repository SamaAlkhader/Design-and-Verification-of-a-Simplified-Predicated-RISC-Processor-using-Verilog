initial begin
    for (int i = 0; i < 1024; i = i+1) instMemory[i] = 0;

    // --- PHASE 1: Forwarding & Priority (I-Type & R-Type) ---
    // Format I: Op(5), Rp(5), Rd(5), Rs(5), Imm(12)
    instMemory[0] = {5'd5, 5'd0, 5'd1, 5'd0, 12'd10};  // ADDI R1, R0, 10
    instMemory[1] = {5'd5, 5'd0, 5'd1, 5'd1, 12'd5};   // ADDI R1, R1, 5  (EX-EX Forward: R1=15)
    instMemory[2] = {5'd5, 5'd0, 5'd1, 5'd1, 12'd2};   // ADDI R1, R1, 2  (Priority Test: R1=17)
    
    // Format R: Op(5), Rp(5), Rd(5), Rs(5), Rt(5), Unused(7)
    instMemory[3] = {5'd0, 5'd0, 5'd2, 5'd1, 5'd1, 7'd0}; // ADD R2, R1, R1 (R2 should be 34)

    // --- PHASE 2: Load-Use & Memory (I-Type) ---
    instMemory[4] = {5'd10, 5'd0, 5'd2, 5'd0, 12'd100}; // SW R2, 100(R0) (Store 34)
    instMemory[5] = {5'd9,  5'd0, 5'd3, 5'd0, 12'd100}; // LW R3, 100(R0) (R3 = 34)
    instMemory[6] = {5'd0,  5'd0, 5'd4, 5'd3, 5'd3, 7'd0}; // ADD R4, R3, R3 (STALL check: R4=68)

    // --- PHASE 3: Predication (Nullify) ---
    instMemory[7] = {5'd5, 5'd0, 5'd5, 5'd0, 12'd0};    // R5 = 0 (Predicate False)
    instMemory[8] = {5'd5, 5'd0, 5'd6, 5'd0, 12'd1};    // R6 = 1 (Predicate True)
    
    // KILLED: Rp=R5 (0), so ADDI R7 shouldn't happen
    instMemory[9] = {5'd5, 5'd5, 5'd7, 5'd0, 12'd100}; 
    // EXECUTED: Rp=R6 (1), so ADDI R8 should happen
    instMemory[10]= {5'd5, 5'd6, 5'd8, 5'd0, 12'd100};  // R8 = 100

    // --- PHASE 4: J-Type Call & Flush ---
    // Format J: Op(5), Rp(5), Offset(22)
    // CALL +3 if R6 != 0. Jumps to PC 14, writes PC+1 (12) to R31.
    instMemory[11] = {5'd12, 5'd6, 22'd3}; 
    
    instMemory[12] = {5'd5, 5'd0, 5'd10, 5'd0, 12'd255}; // Should be FLUSHED
    instMemory[13] = {5'd5, 5'd0, 5'd10, 5'd0, 12'd255}; // Should be FLUSHED
    
    // --- PHASE 5: Indirect Jump (I-Type style) ---
    // Some architectures use I-Type for JRs to keep the Rs field.
    instMemory[14] = {5'd5, 5'd0, 5'd9, 5'd0, 12'd17};  // R9 = 17
    instMemory[15] = {5'd13, 5'd6, 5'd0, 5'd9, 12'd0};  // J Rs (Jump to R9=17)
    instMemory[16] = {5'd5, 5'd0, 5'd10, 5'd0, 12'd255}; // Should be FLUSHED

    instMemory[17] = {5'd5, 5'd0, 5'd11, 5'd31, 12'd0}; // Final Check: R11 = R31 (12)
end
