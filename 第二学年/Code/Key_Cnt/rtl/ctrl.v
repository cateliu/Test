module ctrl(
	Clk,
	Rst_n,
	T,
	TX_Done,
	wrt_en,
	Send_en,
	data
	);
	
	input Clk;
	input Rst_n;
	input [31:0] T;
	input TX_Done;
	input wrt_en;
	
	output reg [7:0] data;
	output reg Send_en;
	reg 	[7:0] T1,T2,T3,T4;
	
	localparam
		WAIT = 4'b0000,
		TX_1 = 4'b0001,
		TX_2 = 4'b0010,
		TX_3 = 4'b0100,
		TX_4 = 4'b1000;
		
	reg	[3:0] state;
	reg 	IsWork;
	
	// 将32位数据分成4份
	always@(posedge Clk or negedge Rst_n)
	if(!Rst_n) begin
		{T1,T2,T3,T4} <= 32'd0;
	end
	else if(wrt_en == 1&& IsWork == 0)begin
		T1 <= T[7:0];
		T2 <= T[15:8];
		T3 <= T[23:16];
		T4 <= T[31:24];
	end

	
	always@(posedge Clk or negedge Rst_n)
	if(!Rst_n)begin
		state 		<= WAIT;
		IsWork 		<= 1'b0;//为0时表示当前正在发送数据
		Send_en 		<= 1'b0;
		data			<=	8'd0;
//		TX_Done		<= 1'b0;
	end
	else begin
		case(state)
			WAIT:
				if(wrt_en == 1)
					begin
						data 			<= T4;		//第一字节数据
						IsWork 		<= 1'b1;	//表示正在发送数据
						Send_en 		<= 1'b1;	//发送使能信号
						state 		<= TX_1;	//状态改变，开始发送第一节数据
					end
				else
					state 			<= WAIT;
			TX_1:
				if(TX_Done === 1)					//如果数据发送完成，TX_Done == 1，那么开始发送第二字节数据
					begin
						state 			<= TX_2;	// 状态改变，开始发送第二节数据
						data 				<= T3;		//第二节数据
						Send_en			<= 1'b1;
					end
				else
					begin
						state 			<= TX_1;	//状态不改变
						Send_en			<= 1'b0;
					end
			TX_2:
				if(TX_Done === 1)					//如果数据发送完成，TX_Done == 1，那么开始发送第三字节数据
					begin
						state 			<= TX_3;	// 状态改变，开始发送第二节数据
						data 				<= T2;		//第二节数据
						Send_en			<= 1'b1;
					end
				else
					begin
						state 			<= TX_2;	//状态不改变
						Send_en			<= 1'b0;
					end
			TX_3:
				if(TX_Done === 1)					//如果数据发送完成，TX_Done == 1，那么开始发送第四字节数据
					begin							
						state 			<= TX_4;	// 状态改变，开始发送第四节数据
						data 				<= T1;	//第四节数据
						Send_en			<= 1'b1;
					end
				else 
					begin
						state 			<= TX_3;	//状态不改变
						Send_en			<= 1'b0;
					end
			TX_4:
				if(TX_Done == 1)					//数据发送完成，TX_Done == 1，四个字节传输完成，状态变为WAIT
					begin
						state 			<= WAIT;	// 状态改变，开始发送第二节数据
						data				<= 8'd0;
						IsWork			<= 1'b0;
						Send_en			<= 1'b0;
					end
				else
				begin
					state 				<= TX_4;	//状态不改变
					Send_en				<= 1'b0;
				end
			default:
				state <= WAIT;
		endcase
	end
	

	
endmodule