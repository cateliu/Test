`timescale 1ns/1ns
`define clk_period 20

module key_time_cnt_tb;
	reg Clk;
	reg Rst_n;
	wire data;
	 
	key_time_cnt time_cnt(
		.Clk			(Clk	),
		.Rst_n		(Rst_n),
		.data_t		(data)
	);
	initial Clk = 1;
	always#(`clk_period/2) Clk=~Clk;
	
	initial begin
		Rst_n = 0;
		#(`clk_period*10) Rst_n = 1;
		
	end


endmodule