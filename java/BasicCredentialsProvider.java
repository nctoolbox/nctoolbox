import org.apache.commons.httpclient.Credentials;
import org.apache.commons.httpclient.UsernamePasswordCredentials;
import org.apache.commons.httpclient.auth.AuthScheme;
import org.apache.commons.httpclient.auth.CredentialsNotAvailableException;
import org.apache.commons.httpclient.auth.CredentialsProvider;

public class BasicCredentialsProvider implements CredentialsProvider {

    /* 
	 * Yeah, only username/password authentication for now.
	 */
    private UsernamePasswordCredentials credentials;

    public BasicCredentialsProvider(String userName, String password) {
        credentials = new UsernamePasswordCredentials(userName,password);
    }

	/*
	 * All we need do to implement the interface is implement
	 * the getCredentials method.
	 */
    public Credentials getCredentials (AuthScheme scheme, String host, int port, boolean proxy) 
    		throws CredentialsNotAvailableException {
		return(credentials);
	}

}