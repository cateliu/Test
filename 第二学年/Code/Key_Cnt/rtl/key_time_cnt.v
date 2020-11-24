module key_time_cnt(
	Clk,
	Rst_n,
	key1,
	data_t
);
	input Clk;
	input Rst_n;
	input key1;
	
	output data_t;

	wire wrt_en;
	wire key_flag,key_state1;
	wire[31:0] T;
	

	key_filter key_filter1(
		.Clk			(Clk		),
		.Rst_n		(Rst_n	),
		.key_in		(key1		),
		.key_flag	(key_flag),
		.key_state	(key_state1)
		);
	wire key_cnt;
	assign key_cnt = key_flag&(~key_state1);
//	reg [31:0] cnt;
//	reg key_cnt;
//	always@(posedge Clk or negedge Rst_n)
//	if(!Rst_n)
//		cnt <= 32'd0;
//	else if(cnt == 32'd100_000)begin
//		cnt <= 32'd0;
//		key_cnt <= 1'b1;
//	end
//	else begin
//		cnt <= cnt +1'b1;
//		key_cnt <= 1'b0;
//	end
	
	
	
	counter counter1(
		.Clk			(Clk		), 
		.Rst_n		(Rst_n	), 
		.Cin			(key_cnt	), 
		.T				(T			), 
		.c_en			(wrt_en), 
		.M				(			)
	);
	
	wire Tx_Done;
	wire Send_en;
	wire [7:0] data;
	ctrl ctrl1(
		.Clk			(Clk		),
		.Rst_n		(Rst_n	),
		.T				(T			),
		.TX_Done		(Tx_Done	),
		.wrt_en		(wrt_en	),
		.Send_en		(Send_en),
		.data			(data		)
	);
	
	uart_byte_tx uart_byte_tx1(
		.Clk			(Clk		),
		.Rst_n		(Rst_n	),
		.data_byte	(data		),
		.send_en		(Send_en	),
		.baud_set	(3'd4		),
		
		.Rs232_Tx	(data_t	),
		.Tx_Done		(Tx_Done	),
		.uart_state	()
	);
		
endmodule