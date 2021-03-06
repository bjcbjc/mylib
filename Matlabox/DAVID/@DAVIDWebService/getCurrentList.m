function xReturn = getCurrentList(obj)
%getCurrentList(obj)
%
%     Input:
%   
%     Output:
%       return = (int)

% Build up the argument lists.
values = { ...
   };
names = { ...
   };
types = { ...
   };

% Create the message, make the call, and convert the response into a variable.
soapMessage = createSoapMessage( ...
    'http://service.session.sample', ...
    'getCurrentList', ...
    values,names,types,'document');
response = callSoapService( ...
    obj.endpoint, ...
    'urn:getCurrentList', ...
    soapMessage);
xReturn = parseSoapResponse(response);
