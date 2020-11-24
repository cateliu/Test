`timescale 1ns/1ns

`define clk_period 20

module Key_Cnt_tb;

	reg Clk;
	reg Rst_n;
	
	reg press0;
	reg press1;
	
	wire Data;
	
	
	Key_Cnt key_cnt_tb(
		.Clk		(Clk		),
		.Rst_n	(Rst_n	),
		.key1		(press0	),
		.key2		(press1	),
		.Rs232_Tx(Data		)
	);
	
	initial Clk = 1;
	
	always#(`clk_period/2) Clk = ~Clk;
	
	initial begin
		Rst_n = 1'b0;
		press0 = 0;
		press1 = 0;
		#(`clk_period*10) Rst_n = 1'b1;
		#(`clk_period*10 + 1);	
		
		press0 = 1;
		#(`clk_period*3)
		press0 = 0;
		
		#80_000_000
		
		press0 = 1;
		#(`clk_period*3)
		press0 = 0;
		
		#80_000_000;
		
		press1 = 1;
		#(`clk_period*3)
		press1 = 0;
		
		#80_000_000;
		
		press1 = 1;
		#(`clk_period*3)
		press1 = 0;
		
		#80_000_000;
		$stop;
		
	end
	
	
endmodule