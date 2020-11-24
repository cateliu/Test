`timescale 1ns/1ns

`define clk_period 20

module	ctr_tb;
	
	reg Clk;
	reg Rst_n;
//	reg [31:0]T;
//	wire TX_Done;
//	reg wrt_en;
//	wire Send_en;
//	wire [7:0]data;
	wire data_t;
	
//	uart_byte_tx uart_byte_tx1(
//		.Clk			(Clk	),
//		.Rst_n		(Rst_n	),
//		.data_byte	(data	),
//		.send_en		(Send_en),
//		.baud_set	(3'd4	),
//		
//		.Rs232_Tx	(data_t	),
//		.Tx_Done		(TX_Done),
//		.uart_state	()
//	);
//	
//	ctrl ctrl1(
//		.Clk			(Clk	),
//		.Rst_n		(Rst_n	),
//		.T				(T		),
//		.TX_Done		(TX_Done),
//		.wrt_en		(wrt_en	),
//		.Send_en		(Send_en),
//		.data			(data	)
//	);
	
	
	ctrl_uart	ctrl_uart1(
		.Clk			(Clk		),
		.Rst_n		(Rst_n	),
//		.T				(T			),
//		.wrt_en		(wrt_en	),		
		.data_t			(data_t	)
	); 
	
	
	initial Clk = 1;
	always#(`clk_period/2) Clk = ~Clk;
	
	initial begin
		Rst_n 	= 	1'b0		;
//		T 			= 	32'h0000	;
////		TX_Done	=	1'b0		;
//		wrt_en	=	1'b0		;
		
		#(`clk_period*2+1) Rst_n = 1'b1;
//		#(`clk_period*2)		
//		T = 32'hABCD;
//		wrt_en = 1;
//		#(`clk_period)
//		wrt_en = 0;
////		#(`clk_period*1000)
		#(`clk_period*50e6) Rst_n = 1'b0;
		$stop;
		
	end
endmodule