module key_ctrl(Clk,Rst_n,key_flag0,key_flag1,key_state0,key_state1,Data_byte,Send_en);

	input Clk;
	input Rst_n;
	input key_flag0,key_flag1;
	
	input key_state0,key_state1;
	
	output [7:0]Data_byte;
	output reg Send_en;
	
	reg [7:0]key_cnt;
	
	always@(posedge Clk or negedge Rst_n)
	if(!Rst_n)begin
		key_cnt <= 7'b000_0000;
		Send_en <= 1'b0;
		end
	else if(key_flag0 && !key_state0)begin
		key_cnt <= key_cnt + 1'b1;
		Send_en <= 1'b1;
	end
	else if(key_flag1 && !key_state1)begin
		key_cnt <= key_cnt - 1'b1;
		Send_en <= 1'b1;
	end
	else begin
		key_cnt <= key_cnt;
		Send_en <= 1'b0;
		end
	
	assign Data_byte = key_cnt;
//	assign Send_en = (key_flag0&(!key_state0))^(key_flag1&(!key_state1));
		
endmodule
