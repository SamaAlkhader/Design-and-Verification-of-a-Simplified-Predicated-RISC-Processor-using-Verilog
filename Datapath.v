`timescale 1ns/1ps

// =============================================================================
// TOP-LEVEL DATAPATH
// =============================================================================
module Datapath(
    input wire clk,
    input wire reset,

    // --- 1. FETCH ---
    output wire [31:0] w_IF_PC, w_IF_Instr,
    output wire [1:0]  w_PCSrc,
    output wire [31:0] w_NextPC,
    
    // --- 2. DECODE ---
    output wire [31:0] w_ID_PC, w_ID_Instr,
    output wire [4:0]  w_Rs_ID, w_Rt_ID, w_Rd_ID,
    output wire        w_JumpOccurred,
    output wire        w_PCWrite, w_IF_ID_Write, w_ID_EX_Flush,

    // --- 3. EXECUTE ---
    output wire [31:0] w_Reg1_EX, w_Reg2_EX, w_Imm_EX, w_PC_EX, w_RetAddr_EX,
    output wire [4:0]  w_Rs_EX, w_Rt_EX, w_Rd_EX, w_Opcode_EX,
    output wire [2:0]  w_ALUOp_EX,
    output wire        w_RegWrite_EX, w_MemRead_EX, w_MemWrite_EX, w_WB_EX, 
    output wire        w_ALUSrc_EX, w_nullify_EX,
    output wire [1:0]  w_WBdata_sel_EX,
    output wire [1:0]  w_ForwardA, w_ForwardB, w_ForwardD,
    output wire [31:0] w_ALU_Comb_EX,

    // --- 4. MEMORY ---
    output wire [31:0] w_ALU_Result_MEM, w_StoreData_MEM, w_RetAddr_MEM,
    output wire [4:0]  w_Rd_MEM,
    output wire        w_RegWrite_MEM, w_MemRead_MEM, w_MemWrite_MEM, w_WB_MEM,
    output wire [1:0]  w_WBdata_sel_MEM,

    // --- 5. WRITE-BACK ---
    output wire [31:0] w_WB_Data_Final,
    output wire [4:0]  w_Rd_WB,
    output wire        w_RegWrite_WB
);

    // --- INTERNAL LOGIC ---
    assign w_PCSrc = w_JumpOccurred ? 2'b10 : 2'b00;
    assign w_Rs_ID = w_ID_Instr[16:12];
    assign w_Rt_ID = w_ID_Instr[11:7];
    assign w_Rd_ID = w_ID_Instr[21:17];

    // --- HAZARD DETECTION UNIT ---
    HazardDetectionUnit HDU(
        .ID_EX_MemRead(w_MemRead_EX),
        .ID_EX_Rd(w_Rd_EX),
        .IF_ID_Rs(w_Rs_ID),
        .IF_ID_Rt(w_Rt_ID),
        .IF_ID_Rp(w_ID_Instr[26:22]), 
        .IF_ID_Rd_Input(w_Rd_ID),
        .IF_ID_Opcode(w_ID_Instr[31:27]),
        .PCWrite(w_PCWrite),
        .IF_ID_Write(w_IF_ID_Write),
        .ID_EX_Flush(w_ID_EX_Flush)
    );
    
    // --- FETCH STAGE ---
    Fetch_Stage IF (
        .clk(clk), .reset(reset),
        .PCWrite(w_PCWrite),
        .PCSrc(w_PCSrc), 
        .Target_Addr(w_NextPC),
        .IF_ID_Write(w_IF_ID_Write),
        .flush(w_JumpOccurred),
        .ID_PC(w_ID_PC), .ID_Instruction(w_ID_Instr),
        .IF_PC(w_IF_PC), .IF_Instr(w_IF_Instr)
    );

    // --- DECODE STAGE ---
    Decode_Stage ID (
        .clk(clk), .reset(reset),
        .ID_EX_Flush(w_ID_EX_Flush),
        .ID_PC(w_ID_PC), .ID_Instruction(w_ID_Instr),
        .RegWr_WB(w_RegWrite_WB), .Rd_WB(w_Rd_WB), .WBbus(w_WB_Data_Final),
        .ALU_EX(w_ALU_Comb_EX), .ALU_MEM(w_ALU_Result_MEM),
        .Rd_EX_in(w_Rd_EX), .Rd_MEM_in(w_Rd_MEM),
        .RegWr_EX_in(w_RegWrite_EX), .RegWr_MEM_in(w_RegWrite_MEM),
        .Reg1_EX(w_Reg1_EX), .Reg2_EX(w_Reg2_EX), .Imm_EX(w_Imm_EX),
        .PC_EX(w_PC_EX), .RetAddr_EX(w_RetAddr_EX), .Rd_EX(w_Rd_EX),
        .Rs_EX(w_Rs_EX), .Rt_EX(w_Rt_EX), 
        .Opcode_EX(w_Opcode_EX), .ALUOp_EX(w_ALUOp_EX), .RegWrite_EX(w_RegWrite_EX),
        .MemRead_EX(w_MemRead_EX), .MemWrite_EX(w_MemWrite_EX), .WB_EX(w_WB_EX),
        .ALUSrc_EX(w_ALUSrc_EX), .WBdata_sel_EX(w_WBdata_sel_EX), .nullify_EX(w_nullify_EX),
        .NextPC(w_NextPC), .JumpOccurred(w_JumpOccurred)
    );

    // --- EXECUTE STAGE ---
    Execute_Stage EX (
        .clk(clk), .reset(reset),
        .Reg1(w_Reg1_EX), .Reg2(w_Reg2_EX), .Imm(w_Imm_EX),
        .Rs_in(w_Rs_EX), .Rt_in(w_Rt_EX), .Rd_in(w_Rd_EX), .RetAddr_in(w_RetAddr_EX),
        .ALUControl(w_ALUOp_EX), .ALUSrc(w_ALUSrc_EX),
        .RegWrite_in(w_RegWrite_EX), .MemRead_in(w_MemRead_EX), .MemWrite_in(w_MemWrite_EX),
        .EX_MEM_ALU_Result(w_ALU_Result_MEM), .MEM_WB_Data_Reg(w_WB_Data_Final),    
        .EX_MEM_Rd(w_Rd_MEM), .MEM_WB_Rd(w_Rd_WB),
        .EX_MEM_RegWrite(w_RegWrite_MEM), .MEM_WB_RegWrite(w_RegWrite_WB),
        
		.WB_in(w_WB_EX),                
        .WBdata_sel(w_WBdata_sel_EX),    	
        .nullify(w_nullify_EX),         	
		
        // Forwarding Selector Connections
        .forwardA(w_ForwardA), .forwardB(w_ForwardB), .forwardD(w_ForwardD),
		
		.RetAddr_MEM(w_RetAddr_MEM),
        .WB_MEM(w_WB_MEM),               
        .WBdata_sel_MEM(w_WBdata_sel_MEM),
        .ALU_Result_Comb(w_ALU_Comb_EX), 
        .ALU_Result_MEM(w_ALU_Result_MEM), .StoreData_MEM(w_StoreData_MEM),
        .Rd_MEM(w_Rd_MEM), .RegWrite_MEM(w_RegWrite_MEM),
        .MemRead_MEM(w_MemRead_MEM), .MemWrite_MEM(w_MemWrite_MEM)
    );

    // --- MEMORY STAGE ---
    Memory_Stage MEM (
        .clk(clk), .reset(reset),
        .ALU_Result(w_ALU_Result_MEM), .StoreData(w_StoreData_MEM),
        .RetAddr_in(w_RetAddr_MEM), .WBdata_sel(w_WBdata_sel_MEM),
		.MemRead(w_MemRead_MEM), .MemWrite(w_MemWrite_MEM),
        .RegWrite_in(w_RegWrite_MEM), .Rd_in(w_Rd_MEM),
        .WB_Data_final(w_WB_Data_Final), .Rd_WB(w_Rd_WB), .RegWrite_WB(w_RegWrite_WB)
    );

endmodule

// =============================================================================
// STAGE MODULES
// =============================================================================

module Fetch_Stage (
    input wire clk, reset,
    input wire PCWrite,
    input wire [1:0] PCSrc,
    input wire [31:0] Target_Addr,
    input wire IF_ID_Write, flush,
	output wire [31:0] IF_PC,
    output wire [31:0] IF_Instr,
    output wire [31:0] ID_PC, ID_Instruction 
	
);

	assign IF_PC = w_PC; 
	assign IF_Instr = w_Instruction;
    wire [31:0] w_PC, w_NextPC, w_Instruction;
    
    PC_Register pc_reg (.clk(clk), .reset(reset),.PCWrite(PCWrite), .next_pc(w_NextPC), .pc(w_PC));
    
    assign w_NextPC = (PCSrc == 2'b10) ? Target_Addr : (w_PC + 32'd1);
    
    InstructionMemory inst_mem (.clk(clk), .addr(w_PC), .instruction(w_Instruction)); 
    
    IF_ID_Register if_id_reg (
        .clk(clk), .reset(reset), 
        .write_en(IF_ID_Write), 
        .flush(flush),
        .PC_in(w_PC), .Instr_in(w_Instruction),
        .PC_out(ID_PC), .Instr_out(ID_Instruction)
    );
endmodule

module Decode_Stage (
    input wire clk, reset,
    input wire ID_EX_Flush,
    input wire [31:0] ID_PC, ID_Instruction,
    input wire RegWr_WB, [4:0] Rd_WB, [31:0] WBbus,
    
    //INPUTS FOR DECODE FORWARDING
    input wire [31:0] ALU_EX, ALU_MEM,
    input wire [4:0]  Rd_EX_in, Rd_MEM_in,
    input wire        RegWr_EX_in, RegWr_MEM_in,

    output wire [31:0] Reg1_EX, Reg2_EX, Imm_EX, PC_EX, RetAddr_EX,
    output wire [4:0]  Rd_EX, Opcode_EX, Rs_EX, Rt_EX,
    output wire [2:0]  ALUOp_EX,
    output wire RegWrite_EX, MemRead_EX, MemWrite_EX, WB_EX, ALUSrc_EX, nullify_EX,
    output wire [1:0]  WBdata_sel_EX,
    output wire [31:0] NextPC,
    output wire JumpOccurred 
);
    wire [31:0] w_RsData, w_RtData, w_RpData, w_ExtImm;
    wire [4:0]  w_Rs, w_Rt, w_Rd, w_Rp, w_Opcode;
    wire        w_ExtOp, w_ALUSrc, w_MemRd, w_MemWr, w_RegWr, w_nullify;
    wire [1:0]  w_WBdata; wire [2:0] w_ALUOp;

    assign w_Opcode = ID_Instruction[31:27];
    assign w_Rp     = ID_Instruction[26:22];
    assign w_Rd     = ID_Instruction[21:17];
    assign w_Rs     = ID_Instruction[16:12];
    assign w_Rt     = ID_Instruction[11:7];
    wire [4:0] w_ReadAddr2;

    // Mux to select between Rt and Rd as the second source register
    assign w_ReadAddr2 = (w_Opcode == 5'd10) ? w_Rd : w_Rt;

    RegisterFile rf (
        .clk(clk), .curr_pc(ID_PC), .Rs1(w_Rs), .Rs2(w_ReadAddr2), .Rp(w_Rp),
        .Rd(Rd_WB), .RegWr(RegWr_WB), .WBbus(WBbus),
        .Bus1(w_RsData), .Bus2(w_RtData), .BusP(w_RpData)
    );

    // --- PREDICATE (Rp) FORWARDING LOGIC ---
    wire [31:0] w_RealRpData;
    assign w_RealRpData = (w_Rp == 5'd0) ? 32'b0 :
                          ((w_Rp == Rd_EX_in)  && RegWr_EX_in)  ? ALU_EX  :
                          ((w_Rp == Rd_MEM_in) && RegWr_MEM_in) ? ALU_MEM :
                          w_RpData;
                        
    MainControlUnit mcu (.Opcode(w_Opcode), .ExtOp(w_ExtOp), .ALUSrc(w_ALUSrc), .ALUOp(w_ALUOp), .MemRd(w_MemRd), .MemWr(w_MemWr), .WBdata(w_WBdata), .RegWr(w_RegWr));
    
    // Pass the CORRECTED RpData to PC Control Unit
    PCControlUnit pcu (
        .PC(ID_PC), .PCPlus1(ID_PC + 32'd1), .Opcode(w_Opcode), 
        .RpIdx(w_Rp), .RpData(w_RealRpData), .RsData(w_RsData), 
        .Offset(ID_Instruction[21:0]), .NextPC(NextPC), .JumpOccurred(JumpOccurred)
    );

    // Use CORRECTED RpData for nullifying the current instruction
    assign w_nullify = (w_Rp != 5'd0) && (w_RealRpData == 32'b0);
    assign w_ExtImm  = (w_ExtOp) ? {{20{ID_Instruction[11]}}, ID_Instruction[11:0]} : {20'b0, ID_Instruction[11:0]};

    wire [4:0] w_Rd_to_EX = (w_Opcode == 5'd12) ? 5'd31 : w_Rd;

    ID_EX_Register id_ex_reg (
        .CLK(clk), .reset(reset),.ID_EX_Flush(ID_EX_Flush), 
        .RegWrite_in(ID_EX_Flush ? 1'b0 : (w_nullify ? 1'b0 : w_RegWr)), 
        .Mem_read_in(ID_EX_Flush ? 1'b0  : (w_nullify ? 1'b0 : w_MemRd)), 
        .Mem_write_in(ID_EX_Flush ? 1'b0 : (w_nullify ? 1'b0 : w_MemWr)), 
        .WB_in (ID_EX_Flush ? 1'b0 : (w_nullify ? 1'b0 : w_RegWr)), 
        .ALUControl_in(w_ALUOp), .ALUSrc_in(w_ALUSrc),
        .WBdata_in(w_WBdata), .Opcode_in(w_Opcode), .nullify_in(w_nullify),.enable(1'b1),
        .Rd_in(w_Rd_to_EX), .Rs_in(w_Rs), .Rt_in(w_Rt), .Reg1_in(w_RsData), .Reg2_in(w_RtData), 
        .Imm_in(w_ExtImm), .PC_in(ID_PC), .RetAddr_in(ID_PC + 32'd1),
        .RegWrite_out(RegWrite_EX), .Mem_read_out(MemRead_EX), .Mem_write_out(MemWrite_EX),
        .WB_out(WB_EX), .ALUControl_out(ALUOp_EX), .ALUSrc_out(ALUSrc_EX), .WBdata_out(WBdata_sel_EX), 
        .Opcode_out(Opcode_EX), .Rd_out(Rd_EX), .Rs_out(Rs_EX), .Rt_out(Rt_EX),
        .Reg1_out(Reg1_EX), .Reg2_out(Reg2_EX), .Imm_out(Imm_EX), .PC_out(PC_EX), .RetAddr_out(RetAddr_EX), .nullify_out(nullify_EX)
    );
endmodule


module Execute_Stage (
    input wire clk, reset,
    input wire [31:0]  Reg1, Reg2, Imm, RetAddr_in,
    input wire [4:0]   Rs_in, Rt_in, Rd_in,
    input wire [2:0]   ALUControl,
    input wire [1:0]   WBdata_sel,
    input wire         ALUSrc, RegWrite_in, MemRead_in, MemWrite_in, WB_in, nullify,
  
    
    // Pipeline Data Inputs
    input wire [31:0]  EX_MEM_ALU_Result, 
    input wire [31:0]  MEM_WB_Data_Reg,      
    
    // Hazard/Forwarding Inputs
    input wire [4:0]   EX_MEM_Rd, 
    input wire [4:0]   MEM_WB_Rd,
    input wire         EX_MEM_RegWrite, 
    input wire         MEM_WB_RegWrite,	
	
	output wire [1:0] forwardA, 
    output wire [1:0] forwardB, 
    output wire [1:0] forwardD,

    output wire [31:0] ALU_Result_MEM, StoreData_MEM, RetAddr_MEM,
    output wire [31:0] ALU_Result_Comb, 
    output wire [4:0]  Rd_MEM,
    output wire        RegWrite_MEM, MemRead_MEM, MemWrite_MEM, WB_MEM,
    output wire [1:0]  WBdata_sel_MEM
);  

    wire [31:0] alu_in_B, alu_out;
    reg  [31:0] fa_mux_out, fb_mux_out, fd_mux_out;

    // --- FORWARDING UNIT INSTANTIATIONS ---
    // Forwarding for Rs
   Forwarding_Unit fwA (
    .ID_EX_Rs(Rs_in),
    .EX_MEM_Rd(EX_MEM_Rd),
    .EX_MEM_RegWrite(EX_MEM_RegWrite),
    .EX_MEM_MemRead(MemRead_MEM), 
    .MEM_WB_Rd(MEM_WB_Rd),
    .MEM_WB_RegWrite(MEM_WB_RegWrite),
    .Forward(forwardA)
);

    
    // Forwarding for Rt 
    Forwarding_Unit fwB (
        .ID_EX_Rs(Rt_in), 
        .EX_MEM_Rd(EX_MEM_Rd), .EX_MEM_RegWrite(EX_MEM_RegWrite), .EX_MEM_MemRead(MemRead_MEM),
        .MEM_WB_Rd(MEM_WB_Rd), .MEM_WB_RegWrite(MEM_WB_RegWrite),
        .Forward(forwardB)
    );

    // Forwarding for Rd (Data to store for SW)
    Forwarding_Unit fwD (
        .ID_EX_Rs(Rd_in), 
        .EX_MEM_Rd(EX_MEM_Rd), .EX_MEM_RegWrite(EX_MEM_RegWrite), .EX_MEM_MemRead(MemRead_MEM),
        .MEM_WB_Rd(MEM_WB_Rd), .MEM_WB_RegWrite(MEM_WB_RegWrite),
        .Forward(forwardD)
    );

    
    always @(*) begin
        case(forwardA)
            2'b01:   fa_mux_out = EX_MEM_ALU_Result; // From EX/MEM
            2'b10:   fa_mux_out = MEM_WB_Data_Reg;   // From MEM/WB
            default: fa_mux_out = Reg1;              // From RegFile
        endcase
    end

    always @(*) begin
        case(forwardB)
            2'b01:   fb_mux_out = EX_MEM_ALU_Result;
            2'b10:   fb_mux_out = MEM_WB_Data_Reg;
            default: fb_mux_out = Reg2;
        endcase
    end

    // Use ALUSrc to decide between Rt_forwarded or the Immediate
    assign alu_in_B = ALUSrc ? Imm : fb_mux_out;
    
    // Store data (for SW) also needs to be forwarded
    always @(*) begin
        case(forwardD)
            2'b01:   fd_mux_out = EX_MEM_ALU_Result;
            2'b10:   fd_mux_out = MEM_WB_Data_Reg;
            default: fd_mux_out = Reg2; 
        endcase
    end

    ALU alu_unit (
        .A(fa_mux_out), 
        .B(alu_in_B), 
        .ALUControl(ALUControl), 
        .ALUResult(alu_out)
    );

    // Assign the combinational ALU output to the module port
    assign ALU_Result_Comb = alu_out;

    // Pipeline Register
    EX_MEM_Register ex_mem_reg (
        .CLK(clk), .reset(reset), 
        .RegWrite_in(RegWrite_in & ~nullify), 
        .MemRead_in(MemRead_in & ~nullify), 
        .MemWrite_in(MemWrite_in & ~nullify), 
        .WB_in(WB_in), 
        .WBdata_in(WBdata_sel),
        .ALU_Result_in(alu_out), 
        .StoreData_in(fd_mux_out),
        .RetAddr_in(RetAddr_in), 
        .Rd_in(Rd_in), 
        .NPC_in(32'b0),
        .RegWrite_out(RegWrite_MEM), .MemRead_out(MemRead_MEM), .MemWrite_out(MemWrite_MEM), .WB_out(WB_MEM), 
        .WBdata_out(WBdata_sel_MEM), .ALU_Result_out(ALU_Result_MEM), .StoreData_out(StoreData_MEM), 
        .RetAddr_out(RetAddr_MEM), .NPC_out(), .Rd_out(Rd_MEM)
    );
endmodule

module Memory_Stage (
    input wire clk, reset,
    input wire [31:0] ALU_Result, StoreData, RetAddr_in,
    input wire [4:0] Rd_in,
    input wire MemRead, MemWrite, RegWrite_in, [1:0] WBdata_sel,
    output wire [31:0] WB_Data_final,
    output wire [4:0] Rd_WB,
    output wire RegWrite_WB
);
    wire [31:0] w_MemData, w_MuxOut;
    DataMemory dmem (.clk(clk), .MemRead(MemRead), .MemWrite(MemWrite), .Address(ALU_Result), .WriteData(StoreData), .ReadData(w_MemData));
    assign w_MuxOut = (WBdata_sel == 2'b01) ? w_MemData : (WBdata_sel == 2'b10) ? RetAddr_in : ALU_Result;
    MEM_WB_Register mem_wb_reg (.clk(clk), .reset(reset), .RegWrite_in(RegWrite_in), .Data_in(w_MuxOut), .Rd_in(Rd_in), .RegWrite_out(RegWrite_WB), .Data_out(WB_Data_final), .Rd_out(Rd_WB));
endmodule

// =============================================================================
// PIPELINE REGISTERS
// =============================================================================

module IF_ID_Register (input clk, reset, write_en, flush, [31:0] PC_in, Instr_in, output reg [31:0] PC_out, Instr_out);
    always @(posedge clk or posedge reset) if (reset || flush) {PC_out, Instr_out} <= 0; else if (write_en) {PC_out, Instr_out} <= {PC_in, Instr_in};
endmodule

module PC_Register (input clk, reset, PCWrite, [31:0] next_pc, output reg [31:0] pc);
    always @(posedge clk or posedge reset)begin 
    if (reset) pc <= 0;
    // only enable update if not stalling
    else if (PCWrite)
    pc <= next_pc;
    end
endmodule

module ID_EX_Register (
    input wire CLK, reset, ID_EX_Flush, RegWrite_in, Mem_read_in, Mem_write_in, WB_in, ALUSrc_in, nullify_in,
    input wire enable, 
    input wire [1:0] WBdata_in, input wire [2:0] ALUControl_in, input wire [4:0] Opcode_in, Rd_in, Rs_in, Rt_in,
    input wire [31:0] Reg1_in, Reg2_in, Imm_in, PC_in, RetAddr_in,
    output reg RegWrite_out, Mem_read_out, Mem_write_out, WB_out, ALUSrc_out, nullify_out,
    output reg [1:0] WBdata_out, output reg [2:0] ALUControl_out, output reg [4:0] Opcode_out, Rd_out, Rs_out, Rt_out,
    output reg [31:0] Reg1_out, Reg2_out, Imm_out, PC_out, RetAddr_out
);
always @(posedge CLK or posedge reset) begin
    if (reset || ID_EX_Flush) begin
        RegWrite_out <= 0; Mem_read_out <= 0; Mem_write_out <= 0; WB_out <= 0; ALUSrc_out <= 0; nullify_out <= 0;
        WBdata_out <= 0; ALUControl_out <= 0; Opcode_out <= 0; Rd_out <= 0; Rs_out <= 0; Rt_out <= 0;
        Reg1_out <= 0; Reg2_out <= 0; Imm_out  <= 0; PC_out   <= 0; RetAddr_out <= 0;
    end
    else if (enable) begin
        RegWrite_out <= RegWrite_in; Mem_read_out <= Mem_read_in; Mem_write_out <= Mem_write_in; WB_out <= WB_in; ALUSrc_out <= ALUSrc_in; nullify_out <= nullify_in;
        WBdata_out <= WBdata_in; ALUControl_out <= ALUControl_in; Opcode_out <= Opcode_in; Rd_out <= Rd_in; Rs_out <= Rs_in; Rt_out <= Rt_in;
        Reg1_out <= Reg1_in; Reg2_out <= Reg2_in; Imm_out  <= Imm_in; PC_out   <= PC_in; RetAddr_out <= RetAddr_in;
    end
end
endmodule

module EX_MEM_Register (
    input wire CLK, reset, RegWrite_in, MemRead_in, MemWrite_in, WB_in, [1:0] WBdata_in,
    input wire [31:0] ALU_Result_in, StoreData_in, RetAddr_in, NPC_in, input wire [4:0] Rd_in,
    output reg RegWrite_out, MemRead_out, MemWrite_out, WB_out, 
    output reg [1:0] WBdata_out,  
    output reg [31:0] ALU_Result_out, StoreData_out, RetAddr_out, NPC_out,
    output reg [4:0] Rd_out      
);
    always @(posedge CLK or posedge reset) begin
        if (reset) begin
            {RegWrite_out, MemRead_out, MemWrite_out, WB_out} <= 0;
            {WBdata_out, ALU_Result_out, StoreData_out, RetAddr_out, NPC_out, Rd_out} <= 0;
        end else begin
            RegWrite_out <= RegWrite_in; MemRead_out <= MemRead_in; MemWrite_out <= MemWrite_in;
            WB_out <= WB_in; WBdata_out <= WBdata_in; ALU_Result_out <= ALU_Result_in;
            StoreData_out <= StoreData_in; RetAddr_out <= RetAddr_in; NPC_out <= NPC_in; Rd_out <= Rd_in;
        end
    end
endmodule

module MEM_WB_Register (
    input wire clk, reset, RegWrite_in, [31:0] Data_in, [4:0] Rd_in, 
    output reg RegWrite_out, 
    output reg [31:0] Data_out, 
    output reg [4:0] Rd_out      
);
    always @(posedge clk or posedge reset) begin
        if (reset) begin
            RegWrite_out <= 0; Data_out <= 0; Rd_out <= 0;
        end else begin
            RegWrite_out <= RegWrite_in; Data_out <= Data_in; Rd_out <= Rd_in;
        end
    end
endmodule

// =============================================================================
// SUB-COMPONENTS
// =============================================================================

module HazardDetectionUnit(
    input wire ID_EX_MemRead,
    input wire [4:0] ID_EX_Rd,
    input wire [4:0] IF_ID_Rs,
    input wire [4:0] IF_ID_Rt,
    input wire [4:0] IF_ID_Rp,
    input wire [4:0] IF_ID_Rd_Input, 
    input wire [4:0] IF_ID_Opcode,
    output reg PCWrite,
    output reg IF_ID_Write,
    output reg ID_EX_Flush
    );
    
    always @(*) begin 
        PCWrite = 1; IF_ID_Write = 1; ID_EX_Flush = 0;
        
        if (ID_EX_MemRead && (ID_EX_Rd != 0)) begin
            if ((ID_EX_Rd == IF_ID_Rs) || (ID_EX_Rd == IF_ID_Rt) || (ID_EX_Rd == IF_ID_Rp)) begin
                PCWrite = 0; IF_ID_Write = 0; ID_EX_Flush = 1;
            end
            else if ((IF_ID_Opcode == 5'd10) && (ID_EX_Rd == IF_ID_Rd_Input)) begin
                PCWrite = 0; IF_ID_Write = 0; ID_EX_Flush = 1;
            end
        end
    end
endmodule
        
module Forwarding_Unit (
    input [4:0] ID_EX_Rs,
    input [4:0] EX_MEM_Rd,
    input EX_MEM_RegWrite,
    input EX_MEM_MemRead,
    input [4:0] MEM_WB_Rd,
    input MEM_WB_RegWrite,
    output reg [1:0] Forward
);  
    always @(*) begin
    Forward = 2'b00;
    
    if (EX_MEM_RegWrite && !EX_MEM_MemRead && (EX_MEM_Rd != 0) && (EX_MEM_Rd == ID_EX_Rs)) 
        Forward = 2'b01;
    else if (MEM_WB_RegWrite && (MEM_WB_Rd != 0) && (MEM_WB_Rd == ID_EX_Rs)) 
        Forward = 2'b10;
    end
endmodule

module ALU (input [31:0] A, B, input [2:0] ALUControl, output reg [31:0] ALUResult);
    always @(*) begin
        case (ALUControl)
            3'b000: ALUResult = A + B;
            3'b001: ALUResult = A - B;
            3'b010: ALUResult = A & B;
            3'b011: ALUResult = A | B;
            3'b100: ALUResult = ~(A | B);
            default: ALUResult = 32'b0;
        endcase
    end
endmodule

module RegisterFile(
    input clk, 
    input [4:0] Rs1, Rs2, Rp, Rd, 
    input [31:0] curr_pc,     
    input RegWr, 
    input [31:0] WBbus, 
    output [31:0] Bus1, Bus2, BusP
);
    reg [31:0] registers [0:31];
    integer i;

    initial begin
        for (i = 0; i < 32; i = i + 1) registers[i] = (i < 17) ? i : 32'h0;
        registers[0] = 0; 
    end
    
    wire bypass_Rs1 = (Rs1 == Rd) && RegWr && (Rs1 != 0) && (Rs1 != 5'd30);
    wire bypass_Rs2 = (Rs2 == Rd) && RegWr && (Rs2 != 0);
    wire bypass_Rp  = (Rp  == Rd) && RegWr && (Rp  != 0);

    assign Bus1 = (Rs1 == 5'd0)  ? 32'b0 : 
                  (Rs1 == 5'd30) ? curr_pc : 
                  (bypass_Rs1)   ? WBbus : 
                  registers[Rs1];

    assign Bus2 = (Rs2 == 5'd0)  ? 32'b0 : 
                  (bypass_Rs2)   ? WBbus : 
                  registers[Rs2];

    assign BusP = (Rp == 5'd0)   ? 32'b0 : 
                  (bypass_Rp)    ? WBbus : 
                  registers[Rp];

    always @(posedge clk) begin
        if (RegWr && (Rd != 5'd0) && (Rd != 5'd30)) begin
            registers[Rd] <= WBbus;
        end
        registers[30] <= curr_pc;
    end
endmodule

module InstructionMemory (input clk, [31:0] addr, output [31:0] instruction);
    reg [31:0] instMemory [0:1023];
    
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
	
    assign instruction = instMemory[addr];
endmodule



module DataMemory (input clk, MemRead, MemWrite, [31:0] Address, WriteData, output reg [31:0] ReadData);
    reg [31:0] memory [0:1023];
    
    initial  begin
        for (integer i = 0; i < 1024; i = i + 1) 
            memory[i] = 32'd0;
        memory[10]= 32'd35; 
        memory[11]= 32'd8;
        memory[20]= 32'd9;
    end
    
    always @(*) begin     
        #1
        if (MemRead) ReadData <= memory[Address[9:0]];
    end
    
    always @(posedge clk) begin
        #1
        if (MemWrite) begin
        memory[Address[9:0]] <= WriteData;
    end
    end
endmodule

module MainControlUnit (input [4:0] Opcode, output reg ExtOp, ALUSrc, MemRd, MemWr, RegWr, output reg [1:0] WBdata, output reg [2:0] ALUOp);
    always @(*) begin
        {ExtOp, ALUSrc, MemRd, MemWr, RegWr} = 5'b0;
        WBdata = 2'b00; ALUOp = 3'b000;
        case(Opcode)
            5'd0:  begin ALUOp=3'b000; RegWr=1; end // ADD
            5'd1:  begin ALUOp=3'b001; RegWr=1; end // SUB
            5'd2:  begin ALUOp=3'b011; RegWr=1; end // OR
            5'd3:  begin ALUOp=3'b100; RegWr=1; end // NOR
            5'd4:  begin ALUOp=3'b010; RegWr=1; end // AND
            5'd5:  begin ALUOp=3'b000; ALUSrc=1; RegWr=1; ExtOp=1; end // ADDI
            5'd6:  begin ALUOp=3'b011; ALUSrc=1; RegWr=1; end // ORI
            5'd7:  begin ALUOp=3'b100; ALUSrc=1; RegWr=1; end // NORI
            5'd8:  begin ALUOp=3'b010; ALUSrc=1; RegWr=1; end // ANDI
            5'd9:  begin ALUOp=3'b000; ALUSrc=1; MemRd=1; WBdata=2'b01; RegWr=1; ExtOp=1; end // LW
            5'd10: begin ALUOp=3'b000; ALUSrc=1; MemWr=1; ExtOp=1; end // SW
            5'd12: begin RegWr=1; WBdata=2'b10; end // CALL
            default: ;
        endcase
    end
endmodule

module PCControlUnit(
    input [31:0] PC, PCPlus1, RpData, RsData, 
    input [4:0] Opcode, RpIdx,
    input [21:0] Offset, 
    output reg [31:0] NextPC, 
    output reg JumpOccurred
);
    wire [31:0] SignExtOffset = {{10{Offset[21]}}, Offset};
    wire predicate_met = (RpIdx == 5'd0) || (RpData != 32'd0);

    always @(*) begin
        NextPC = PCPlus1; 
        JumpOccurred = 0;
        if (predicate_met) begin
            case (Opcode)
                5'd11: begin NextPC = PC + SignExtOffset; JumpOccurred = 1; end 
                5'd12: begin NextPC = PC + SignExtOffset; JumpOccurred = 1; end 
                5'd13: begin NextPC = RsData;             JumpOccurred = 1; end
                default: ;
            endcase
        end
    end
endmodule

module tb_Datapath();
    

    // Iterators for simulation logic
    integer i, reg_idx;

    // --- 1. CLOCK & RESET SIGNALS ---
    reg clk;
    reg reset;

    // --- 2. FETCH STAGE WIRES ---
    wire [31:0] w_IF_PC, w_IF_Instr;
    wire [1:0]  w_PCSrc;
    wire [31:0] w_NextPC;
    
    // --- 3. DECODE STAGE WIRES ---
    wire [31:0] w_ID_PC, w_ID_Instr;
    wire [4:0]  w_Rs_ID, w_Rt_ID, w_Rd_ID;
    wire        w_JumpOccurred;
    wire        w_PCWrite, w_IF_ID_Write, w_ID_EX_Flush;

    // --- 4. EXECUTE STAGE WIRES ---
    wire [31:0] w_Reg1_EX, w_Reg2_EX, w_Imm_EX, w_PC_EX, w_RetAddr_EX;
    wire [4:0]  w_Rs_EX, w_Rt_EX, w_Rd_EX, w_Opcode_EX;
    wire [2:0]  w_ALUOp_EX;
    wire        w_RegWrite_EX, w_MemRead_EX, w_MemWrite_EX, w_WB_EX;
    wire        w_ALUSrc_EX, w_nullify_EX;
    wire [1:0]  w_WBdata_sel_EX;
    wire [1:0]  w_ForwardA, w_ForwardB, w_ForwardD;
    wire [31:0] w_ALU_Comb_EX;

    // --- 5. MEMORY STAGE WIRES ---
    wire [31:0] w_ALU_Result_MEM, w_StoreData_MEM, w_RetAddr_MEM;
    wire [4:0]  w_Rd_MEM;
    wire        w_RegWrite_MEM, w_MemRead_MEM, w_MemWrite_MEM, w_WB_MEM;
    wire [1:0]  w_WBdata_sel_MEM;

    // --- 6. WRITE-BACK STAGE WIRES ---
    wire [31:0] w_WB_Data_Final;
    wire [4:0]  w_Rd_WB;
    wire        w_RegWrite_WB;

    // --- 7. DATAPATH INSTANTIATION (UUT) ---
    Datapath uut (
        .clk(clk),
        .reset(reset),

        // Fetch
        .w_IF_PC(w_IF_PC),
        .w_IF_Instr(w_IF_Instr),
        .w_PCSrc(w_PCSrc),
        .w_NextPC(w_NextPC),

        // Decode
        .w_ID_PC(w_ID_PC),
        .w_ID_Instr(w_ID_Instr),
        .w_Rs_ID(w_Rs_ID),
        .w_Rt_ID(w_Rt_ID),
        .w_Rd_ID(w_Rd_ID),
        .w_JumpOccurred(w_JumpOccurred),
        .w_PCWrite(w_PCWrite),
        .w_IF_ID_Write(w_IF_ID_Write),
        .w_ID_EX_Flush(w_ID_EX_Flush),

        // Execute
        .w_Reg1_EX(w_Reg1_EX),
        .w_Reg2_EX(w_Reg2_EX),
        .w_Imm_EX(w_Imm_EX),
        .w_PC_EX(w_PC_EX),
        .w_RetAddr_EX(w_RetAddr_EX),
        .w_Rs_EX(w_Rs_EX),
        .w_Rt_EX(w_Rt_EX),
        .w_Rd_EX(w_Rd_EX),
        .w_Opcode_EX(w_Opcode_EX),
        .w_ALUOp_EX(w_ALUOp_EX),
        .w_RegWrite_EX(w_RegWrite_EX),
        .w_MemRead_EX(w_MemRead_EX),
        .w_MemWrite_EX(w_MemWrite_EX),
        .w_WB_EX(w_WB_EX),
        .w_ALUSrc_EX(w_ALUSrc_EX),
        .w_nullify_EX(w_nullify_EX),
        .w_WBdata_sel_EX(w_WBdata_sel_EX),
        .w_ForwardA(w_ForwardA),
        .w_ForwardB(w_ForwardB),
        .w_ForwardD(w_ForwardD),
        .w_ALU_Comb_EX(w_ALU_Comb_EX),

        // Memory
        .w_ALU_Result_MEM(w_ALU_Result_MEM),
        .w_StoreData_MEM(w_StoreData_MEM),
        .w_RetAddr_MEM(w_RetAddr_MEM),
        .w_Rd_MEM(w_Rd_MEM),
        .w_RegWrite_MEM(w_RegWrite_MEM),
        .w_MemRead_MEM(w_MemRead_MEM),
        .w_MemWrite_MEM(w_MemWrite_MEM),
        .w_WB_MEM(w_WB_MEM),
        .w_WBdata_sel_MEM(w_WBdata_sel_MEM),

        // Write-Back
        .w_WB_Data_Final(w_WB_Data_Final),
        .w_Rd_WB(w_Rd_WB),
        .w_RegWrite_WB(w_RegWrite_WB)
    );

    always #5 clk = ~clk;

    initial begin
        clk = 0; reset = 1;
        #20 reset = 0;
        
        $display("\n=================================================================================");
        $display("                    PIPELINED PROCESSOR DETAILED TRACE                           ");
        $display("=================================================================================");

        for (i = 0; i < 29; i = i + 1) begin
             @(posedge clk);
             #1; // Wait for logic to settle     
             
             $display("\n[CYCLE %2d] -------------------------------------------------------------", i);
             
             // --- FETCH STAGE ---
             $display("  [IF STAGE]");
             $display("    PC Current : %h", uut.IF.w_PC);
             $display("    Next PC    : %h", uut.IF.w_NextPC);
             $display("    Instr Raw  : %h", uut.IF.w_Instruction);
             $display("    Stall/Flush: PCWrite=%b | IF/ID_Write=%b | Flush=%b", 
                      uut.w_PCWrite, uut.w_IF_ID_Write, uut.w_JumpOccurred);

             // --- DECODE STAGE ---
             $display("  [ID STAGE] PC: %h", uut.ID.ID_PC);
             $display("    Decoding   : Op=%2d | Rs=%2d | Rt=%2d | Rd=%2d | Rp=%2d", 
                      uut.ID.w_Opcode, uut.ID.w_Rs, uut.ID.w_Rt, uut.ID.w_Rd, uut.ID.w_Rp);
             $display("    Values Read: Reg[%2d]=%h | Reg[%2d]=%h | Reg[%2d]=%h", 
                      uut.ID.w_Rs, uut.ID.w_RsData, uut.ID.w_ReadAddr2, uut.ID.w_RtData, uut.ID.w_Rp, uut.ID.w_RpData);
             $display("    Control    : Predicate_Kill=%b | JumpOccurred=%b | ID_EX_Flush=%b", 
                      uut.ID.w_nullify, uut.w_JumpOccurred, uut.w_ID_EX_Flush);

             // --- EXECUTE STAGE ---
             $display("  [EX STAGE]");
             $display("    ALU Ops    : OpCode=%3b | ALUSrc=%b", uut.EX.ALUControl, uut.EX.ALUSrc);
             $display("    Forwarding : FwdA=%2b (Val=%h) | FwdB=%2b (Val=%h) | FwdD=%2b (Val=%h)", 
                      uut.EX.forwardA, uut.EX.fa_mux_out, 
                      uut.EX.forwardB, uut.EX.fb_mux_out,
                      uut.EX.forwardD, uut.EX.fd_mux_out);
             $display("    ALU Calc   : A=%h | B=%h | Result=%h", 
                      uut.EX.fa_mux_out, uut.EX.alu_in_B, uut.EX.alu_out);
             $display("    Output     : Dest Reg=R%2d | WriteEn=%b | MemRead=%b", 
                      uut.EX.Rd_in, uut.EX.RegWrite_in, uut.EX.MemRead_in);

             // --- MEMORY STAGE ---
             $display("  [MEM STAGE]");
             $display("    Access     : Addr=%h | WriteData=%h", uut.MEM.ALU_Result, uut.MEM.StoreData);
             $display("    Control    : MemRead=%b | MemWrite=%b", uut.MEM.MemRead, uut.MEM.MemWrite);
             $display("    Read Result: %h (from memory)", uut.MEM.dmem.ReadData);

             // --- WRITE BACK STAGE ---
             $display("  [WB STAGE]");
             $display("    WriteBack  : Writing Data %h to Register R%2d", uut.w_WB_Data_Final, uut.w_Rd_WB);
             $display("    Enable     : %b", uut.w_RegWrite_WB);
        end

        // PRINT ALL REGISTERS
        $display("\n================ FINAL REGISTER FILE STATUS ================");
        $display("Reg |  Value     |  Reg |  Value     |  Reg |  Value     ");
        $display("------------------------------------------------------------");
        for (reg_idx = 0; reg_idx < 32; reg_idx = reg_idx + 4) begin
            $display("R%2d : %h | R%2d : %h | R%2d : %h | R%2d : %h", 
                reg_idx,   uut.ID.rf.registers[reg_idx],
                reg_idx+1, uut.ID.rf.registers[reg_idx+1],
                reg_idx+2, uut.ID.rf.registers[reg_idx+2],
                reg_idx+3, uut.ID.rf.registers[reg_idx+3]);
        end
        $display("============================================================\n");
        // --- FINAL DATA MEMORY STATUS ---
        $display("\n================ FINAL DATA MEMORY STATUS ================");
        $display("Address | Value      | Address | Value");
        $display("----------------------------------------------------------");
        $display("Addr 10 : %h (%d)", uut.MEM.dmem.memory[10], uut.MEM.dmem.memory[10]);
        $display("Addr 20 : %h (%d)", uut.MEM.dmem.memory[20], uut.MEM.dmem.memory[20]);
        $display("==========================================================\n");
        $finish;
    end
endmodule
